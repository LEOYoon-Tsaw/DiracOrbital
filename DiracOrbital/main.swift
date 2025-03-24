//
//  main.swift
//  DiracOrbital
//
//  Created by Leo Liu on 3/23/25.
//

import Foundation
import Numerics

private struct IntermediateResults {
    var zalpha: Double?
    var gamma: Double?
    var lambda: Double?
    var gfLambda: Double?
    var relativeEnergy: Double?
    var normalization: Double?
}

private func sgn(_ x: Int) -> Double {
    return if x >= 0 { 1.0 } else { -1.0 }
}

// --- Helper: Double Factorial ---
// Computes n!! for odd n (e.g., (2m-1)!!)
private func doubleFactorial(_ n: Int) -> Double {
    if n <= 1 {
        return 1.0
    } else {
        return Double(n) * doubleFactorial(n - 2)
    }
}

// --- Helper: Factorial ---
// Computes n! for odd n (e.g., (2m-1)!)
private func factorial(_ n: Int) -> Double {
    if n <= 1 {
        return 1.0
    } else {
        return Double(n) * factorial(n - 1)
    }
}

// --- Associated Legendre Polynomials ---
// Uses recurrence relations:
// P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)
// P_{m+1}^m(x) = x (2m+1) P_m^m(x)
// P_l^m(x) = ((2l-1)x P_{l-1}^m(x) - (l+m-1) P_{l-2}^m(x))/(l-m)  for l > m
private func associatedLegendre(l: Int, m: Int, x: Double) -> Double {
    let absM = abs(m)
    let sign = absM % 2 == 0 ? 1.0 : -1.0
    if l < absM {
        return 0.0
    }
    if m < 0 {
        return sign * factorial(l - absM) / factorial(l + absM) * associatedLegendre(l: l, m: absM, x: x)
    }
    
    // Base: P_absM^absM(x)
    var p_mm = 1.0
    if absM > 0 {
        let factor = 1.0 - x * x
        p_mm = sign * doubleFactorial(2 * absM - 1) * pow(factor, Double(absM) / 2)
    }
    
    if l == absM {
        return p_mm
    }
    
    // Next: P_{absM+1}^absM(x)
    let p_mmp1 = x * Double(2 * absM + 1) * p_mm
    if l == absM + 1 {
        return p_mmp1
    }
    
    // Recurrence for l > absM+1
    return ((Double(2 * l - 1) * x * associatedLegendre(l: l-1, m: absM, x: x)) - (Double(l + absM - 1) * associatedLegendre(l: l-2, m: absM, x: x))) / Double(l - absM)
}

// --- Spherical Harmonic Function ---
// Computes Y_l^m(theta, phi) = norm * P_l^m(cos(theta)) * exp(i * m * phi)
private func sphericalHarmonic(l: Int, m: Int, theta: Double, phi: Double) -> Complex<Double> {
    if l < 0 {
        return sphericalHarmonic(l: -l - 1, m: m, theta: theta, phi: phi)
    }
    let absM = abs(m)
    if l > 0 && absM > l {
        return .zero
    }
    // Normalization factor: sqrt((2l+1)/(4π) * ((l-|m|)!/(l+|m|)!))
    let normFactor = sqrt((Double(2 * l + 1) / (4 * Double.pi)) *
                          factorial(l - absM) / factorial(l + absM))
    
    // Associated Legendre polynomial evaluated at cos(theta)
    let legendreValue = associatedLegendre(l: l, m: m, x: cos(theta))
    
    // e^(i m phi)
    let expFactor = Complex.exp(Complex(imaginary: Double(m) * phi))
    
    // Combine all parts: norm * P_l^m(cos(theta)) * exp(i m phi)
    return Complex<Double>(normFactor) * Complex(legendreValue, 0) * expFactor
}

// --- Generalized Laguerre Polynomial ---
// Implements Lₙ^(α)(x) via a recurrence relation.
// (For n = 0: 1, for n = 1: α+1−x, and for n > 1 the recurrence is used.)
private func generalizedLaguerre(n: Int, alpha: Double, x: Double) -> Double {
    if n < 0 {
        return 0.0
    } else if n == 0 {
        return 1.0
    } else if n == 1 {
        return alpha + 1.0 - x
    } else {
       return ((Double(2 * n - 1) + alpha - x) * generalizedLaguerre(n: n-1, alpha: alpha, x: x) - (Double(n - 1) + alpha) * generalizedLaguerre(n: n-2, alpha: alpha, x: x)) / Double(n)
    }
}

private func spinorCoef(a: Int, b: Int) -> Double {
    return sqrt((0.5 * Double(b + 1) + Double(a)) / Double(2 * a + 1))
}

class HydrogenOrbital {
    static let alpha = 0.0072973525643                 // fine structure constant
    static let c = 299792458.0                         // speed of light (m/s)
    static let hbar = 1.054571817e-34                  // reduced Planck constant (J·s)
    static let me = 9.10938356e-31                     // electron mass in kg
    static let rBohr = hbar / (me * c * alpha)         // Bohr radius in m
    static let halfPeriod = Double.pi * rBohr * alpha / c    // half period of global phase change
    static let mechbar = me * c / hbar                 // me * c / hbar
    
    private var intermediateResults = IntermediateResults()
    
    var z: Int
    var n: Int
    var kappa: Int
    var mx2: Int
    
    // Z: nuclear charge number, 1-137
    // n: principal quantum number, >=1
    // kappa: Dirac quantum number, -n...-1 or 1...(n-1)
    // mx2: 2 * m, magnetic quantum number, -2|kappa|+1 to 2|kappa|-1 with stride of 2
    init(z: Int, n: Int, kappa: Int, mx2: Int) {
        guard z >= 1 && z <= 137 else { fatalError("Z must be between 1 and 137") }
        guard n >= 1 else { fatalError("n must be positive") }
        guard (kappa >= -n && kappa <= -1) || (kappa >= 1 && kappa <= n-1) else { fatalError("kappa must be in -n...-1 or 1...(n-1)") }
        guard (mx2 + 1).isMultiple(of: 2) && mx2 > -2 * n && mx2 < 2 * n else { fatalError("mx2 must be an odd integer between -2|kappa|+1 and 2|kappa|-1") }
        self.z = z
        self.n = n
        self.kappa = kappa
        self.mx2 = mx2
    }
    
    private var zalpha: Double {
        if let zalpha = intermediateResults.zalpha {
            return zalpha
        } else {
            let zalpha = Double(self.z) * Self.alpha
            intermediateResults.zalpha = zalpha
            return zalpha
        }
    }
    
    private var gamma: Double {
        if let gamma = intermediateResults.gamma {
            return gamma
        } else {
            let gamma = sqrt(Double(kappa * kappa) - zalpha * zalpha)
            intermediateResults.gamma = gamma
            return gamma
        }
    }
    
    private var nk: Int {
        n - abs(kappa)
    }
    
    private var relativeEnergy: Double {
        if let relativeEnergy = intermediateResults.relativeEnergy {
            return relativeEnergy
        } else {
            let inner = zalpha / (Double(nk) + gamma)
            let relativeEnergy = pow(1 + inner * inner, -0.5)
            intermediateResults.relativeEnergy = relativeEnergy
            return relativeEnergy
        }
    }
    
    private var lambda: Double {
        if let lambda = intermediateResults.lambda {
            return lambda
        } else {
            let lambda = Self.mechbar * sqrt(1 - relativeEnergy * relativeEnergy)
            intermediateResults.lambda = lambda
            return lambda
        }
    }
    
    private var gfLambda: Double {
        if let gfLambda = intermediateResults.gfLambda {
            return gfLambda
        } else {
            let gfLambda = (gamma - Double(kappa) * relativeEnergy) / sqrt(1 - relativeEnergy * relativeEnergy)
            intermediateResults.gfLambda = gfLambda
            return gfLambda
        }
    }
    
    private var normalization: Double {
        if let normalization = intermediateResults.normalization {
            return normalization
        } else {
            let commonFactor = Double(2 * kappa) * (Double(kappa) - gamma)
            var normalization: Double
            if nk != 0 {
                let factorateRatio = Double.gamma(Double(nk)) / Double.gamma(Double(nk + 1) + 2 * gamma)
                let kGammaTerm = Double(kappa) * relativeEnergy / gamma
                normalization = lambda / (Double(nk) + gamma) * factorateRatio * 0.5 * (kGammaTerm * kGammaTerm + kGammaTerm)
            } else {
                normalization = lambda * 0.5 * pow(zalpha / gamma, 2.0) / Double(kappa * kappa) / Double.gamma(1 + 2 * gamma)
            }
            normalization = sqrt(normalization / commonFactor)
            intermediateResults.normalization = normalization
            return normalization
        }
    }
    
    // --- Radial functions ---
    // The radial part is split into two functions g(ρ) and f(ρ) as in the note.
    private func radial(rho: Double) -> (Double, Double) {
        let commonFactor = Double.pow(rho, gamma) * exp(-rho / 2)
        let firstTerm = rho * generalizedLaguerre(n: self.nk - 1, alpha: 2 * self.gamma + 1, x: rho)
        let lastTerm = self.gfLambda * generalizedLaguerre(n: self.nk, alpha: 2 * self.gamma - 1, x: rho)
        let gammaKappa = self.gamma - Double(self.kappa)
        let gFunc = commonFactor * (zalpha * firstTerm + gammaKappa * lastTerm)
        let fFunc = commonFactor * (zalpha * lastTerm + gammaKappa * firstTerm)
        
        return (gFunc, fFunc)
    }
    
    private func spherical(theta: Double, phi: Double) -> (Complex<Double>, Complex<Double>, Complex<Double>, Complex<Double>) {
        // Angular parts: compute the spinor spherical harmonics.
        // The four components use indices as in the note:
        //  • Component 1: Ω(θ,ϕ)_κ^(m – ½)
        //  • Component 2: Ω(θ,ϕ)_κ^(m + ½)
        //  • Component 3: Ω(θ,ϕ)_₋κ^(m – ½)
        //  • Component 4: Ω(θ,ϕ)_₋κ^(m + ½)
        let Y1 = Complex(spinorCoef(a: kappa, b: -mx2)) * sphericalHarmonic(l: kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        let Y2 = Complex(sgn(kappa) * spinorCoef(a: kappa, b: mx2)) * sphericalHarmonic(l: kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        let Y3 = Complex(spinorCoef(a: -kappa, b: -mx2)) * sphericalHarmonic(l: -kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        let Y4 = Complex(sgn(kappa) * spinorCoef(a: -kappa, b: mx2)) * sphericalHarmonic(l: -kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        return (Y1, Y2, Y3, Y4)
    }
    
    func waveFunction(t: Double, r: Double, theta: Double, phi: Double) -> [Complex<Double>] {
        // Radial parts: g(ρ) and f(ρ)
        let rho = 2 * lambda * r
        let (g, f) = radial(rho: rho)
        
        // Time-dependent phase factor exp(–i ε t/ħ).
        let energy = Self.mechbar * Self.c * relativeEnergy
        let globalPhase = Complex<Double>.exp(Complex(imaginary: -energy * t))
        
        // Overall multiplicative factor.
        let commonFactor = globalPhase * Complex(normalization / r)  // (This factor includes normalization, time evolution and r^-1.)
        
        // Spherical parts: Y_κ^m
        let (YPosUp, YPosDn, YNegUp, YNegDn) = spherical(theta: theta, phi: phi)
        // Assemble the four components of the Dirac spinor:
        //   Component 1: g(ρ) · Ω₁
        //   Component 2: – g(ρ) · Ω₂
        //   Component 3: i f(ρ) · Ω₃
        //   Component 4: i f(ρ) · Ω₄
        let posUp = commonFactor * Complex(g) * YPosUp
        let posDn = commonFactor * Complex(-g) * YPosDn
        let negUp = commonFactor * Complex(0, f) * YNegUp
        let negDn = commonFactor * Complex(0, f) * YNegDn
        
        return [posUp, posDn, negUp, negDn]
    }
    
    func waveFunction(t: Double, x: Double, y: Double, z: Double) -> [Complex<Double>] {
        let rxy = Double.hypot(x, y)
        let r = Double.hypot(rxy, z)
        let theta = Double.atan2(y: rxy, x: z)
        let phi = Double.atan2(y: y, x: x)
        return waveFunction(t: t, r: r, theta: theta, phi: phi)
    }
}

// --- Example usage ---
// (These example parameters are illustrative. In practice, one must choose t, r, θ, ϕ
//  as well as the quantum numbers appropriately.)

let startTime = Date()
let orbital = HydrogenOrbital(z: 79, n: 4, kappa: -1, mx2: 1)

var totalProbability = 0.0
for x in (-100)...100 {
    for y in (-100)...100 {
        for z in (-100)...100 {
            let wave = orbital.waveFunction(t: 0, x: Double(x) / 100 * HydrogenOrbital.rBohr, y: Double(y) / 100 * HydrogenOrbital.rBohr, z: Double(z) / 100 * HydrogenOrbital.rBohr)
            var probability = 0.0
            for component in wave {
                probability += (component.conjugate * component).real
            }
            if !probability.isNaN {
                totalProbability += probability * pow(HydrogenOrbital.rBohr / 100, 3)
            }
        }
    }
}
let endTime = Date()

print("Total probability in space is \(totalProbability), calculated in \(endTime.timeIntervalSince(startTime)) seconds.")
