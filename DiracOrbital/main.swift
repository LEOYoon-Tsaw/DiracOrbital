//
//  main.swift
//  DiracOrbital
//
//  Created by Leo Liu on 3/23/25.
//

import Foundation
import Numerics

extension Formatter {
    static let scientific: NumberFormatter = {
        let formatter = NumberFormatter()
        formatter.numberStyle = .scientific
        formatter.positiveFormat = "0.###E+0"
        formatter.exponentSymbol = "e"
        return formatter
    }()
}

private func phaseInPi(_ z: Complex<Double>) -> Double {
    if z.isZero {
        return 0.0
    }
    var phase = z.phase / Double.pi
    if phase < 0 {
        phase += 2.0
    }
    return phase
}

struct BiSpinor: CustomStringConvertible {
    let posiUp: Complex<Double>
    let posiDown: Complex<Double>
    let negaUp: Complex<Double>
    let negaDown: Complex<Double>
    
    var description: String {
        var desp = "ψ1: \(String(format: "%.3e", posiUp.length))∠\(String(format: "%.2f", phaseInPi(posiUp)))π\n"
        desp += "ψ2: \(String(format: "%.3e", posiDown.length))∠\(String(format: "%.2f", phaseInPi(posiDown)))π\n"
        desp += "ψ3: \(String(format: "%.3e", negaUp.length))∠\(String(format: "%.2f", phaseInPi(negaUp)))π\n"
        desp += "ψ4: \(String(format: "%.3e", negaDown.length))∠\(String(format: "%.2f", phaseInPi(negaDown)))π"
        return desp
    }
    
    var density: Double {
        return (posiUp.conjugate * posiUp).real + (posiDown.conjugate * posiDown).real + (negaUp.conjugate * negaUp).real + (negaDown.conjugate * negaDown).real
    }
}

private struct IntermediateResults {
    var zalpha: Double?
    var gamma: Double?
    var lambda: Double?
    var gfLambda: Double?
    var relativeEnergy: Double?
    var energy: Complex<Double>?
    var normalization: Double?
    var spinorCoef: [String: (Double, Int, Int)] = [:]
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

//starting from left to right, (a_0 * x^0 + a_1 * x^1 + ...)
// (For n = 0: 1, for n = 1: α+1−x, and for n > 1 the recurrence is used.)
private func generalizedLaguerre(n: Int, a: Double, x: Double) -> Double {
    if n < 0 {
        return 0
    } else if n == 0 {
        return 1
    } else if n == 1 {
        return a + 1 - x
    } else {
        return ((Double(2 * n - 1) + a - x) * generalizedLaguerre(n: n-1, a: a, x: x) - (Double(n - 1) + a) * generalizedLaguerre(n: n-2, a: a, x: x)) / Double(n)
    }
}

// starting from left to right, (a_0 * x^0 + a_1 * x^1 + ...) * sqrt(1 - x^2)^m
// Uses recurrence relations:
// P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)
// P_{m+1}^m(x) = x (2m+1) P_m^m(x)
// P_l^m(x) = ((2l-1)x P_{l-1}^m(x) - (l+m-1) P_{l-2}^m(x))/(l-m)  for l > m
private func associatedLegendre(l: Int, m: Int, cosTheta: Double, sinTheta: Double) -> Double {
    guard l >= 0 else {
        return associatedLegendre(l: -l-1, m: m, cosTheta: cosTheta, sinTheta: sinTheta)
    }
    guard abs(m) <= l else {
        return 0
    }
    
    let sign = m.isMultiple(of: 2) ? 1.0 : -1.0
    guard m >= 0 else {
        return sign * factorial(l + m) / factorial(l - m) * associatedLegendre(l: l, m: -m, cosTheta: cosTheta, sinTheta: sinTheta)
    }
    
    let ppm = sign * doubleFactorial(2 * m - 1) * Double.pow(sinTheta, Double(m))
    if l == m {
        return ppm
    } else if l == m + 1 {
        return ppm * Double(2 * m + 1) * cosTheta
    } else {
        return (cosTheta * Double(2 * l - 1) * associatedLegendre(l: l-1, m: m, cosTheta: cosTheta, sinTheta: sinTheta) - Double(l + m - 1) * associatedLegendre(l: l-2, m: m, cosTheta: cosTheta, sinTheta: sinTheta)) / Double(l - m)
    }
}

private func sphericalHarmonicNormalization(l: Int, m: Int) -> Double {
    if l >= 0 {
        return sqrt(Double(2 * l + 1) / (4 * Double.pi) * factorial(l - m) / factorial(l + m))
    } else {
        return sqrt(Double(-2 * l - 1) / (4 * Double.pi) * factorial(-l - 1 - m) / factorial(-l - 1 + m))
    }
}

private func spinorCoef(a: Int, b: Int, downSpinor: Bool, positive: Bool) -> (Double, Int, Int) {
    let m: Int
    let ym: Int
    let k: Int
    let sign: Double
    switch (positive, downSpinor) {
    case (true, false):
        m = -b
        ym = (b-1) / 2
        k = a
        sign = 1.0
    case (true, true):
        m = b
        ym = (b+1) / 2
        k = a
        sign = sgn(a)
    case (false, false):
        m = -b
        ym = (b-1) / 2
        k = -a
        sign = 1.0
    case (false, true):
        m = b
        ym = (b+1) / 2
        k = -a
        sign = sgn(a)
    }
    
    var coef = sphericalHarmonicNormalization(l: k, m: ym)
    coef *= sqrt((0.5 * Double(m + 1) + Double(k)) / Double(2 * k + 1))
    coef *= sign
    return (coef, k, ym)
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
    
    private var energy: Complex<Double> {
        if let energy = intermediateResults.energy {
            return energy
        } else {
            let energy = Complex(imaginary: -Self.mechbar * Self.c * relativeEnergy)
            intermediateResults.energy = energy
            return energy
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
            normalization = 2 * lambda * sqrt(normalization / commonFactor)
            intermediateResults.normalization = normalization
            return normalization
        }
    }
    
    // --- Radial functions ---
    // The radial part is split into two functions g(ρ) and f(ρ) as in the note.
    private func radial(rho: Double) -> (Double, Double) {
        let commonFactor = Double.pow(rho, gamma - 1) * exp(-rho / 2)
        let firstTerm = rho * generalizedLaguerre(n: nk - 1, a: 2 * gamma + 1, x: rho)
        let lastTerm = self.gfLambda * generalizedLaguerre(n: nk, a: 2 * gamma - 1, x: rho)
        let gammaKappa = self.gamma - Double(self.kappa)
        let gFunc = commonFactor * (zalpha * firstTerm + gammaKappa * lastTerm)
        let fFunc = commonFactor * (zalpha * lastTerm + gammaKappa * firstTerm)
        return (gFunc, fFunc)
    }
    
    private func sphericalHarmonics(downSpinor: Bool, positive: Bool, thetaPhase: Complex<Double>, phi: Double) -> Complex<Double> {
        let normalization: Double
        let k: Int
        let m: Int
        if let spinorCoefs = intermediateResults.spinorCoef["\(positive)_\(downSpinor)"] {
            (normalization, k, m) = spinorCoefs
        } else {
            (normalization, k, m) = spinorCoef(a: kappa, b: mx2, downSpinor: downSpinor, positive: positive)
            intermediateResults.spinorCoef["\(positive)_\(downSpinor)"] = (normalization, k, m)
        }
        let polynomial = associatedLegendre(l: k, m: m, cosTheta: thetaPhase.real, sinTheta: thetaPhase.imaginary)
        return Complex(normalization * polynomial) * Complex.exp(Complex(imaginary: phi * Double(m)))
    }
    
    private func spherical(theta: Double, phi: Double) -> (Complex<Double>, Complex<Double>, Complex<Double>, Complex<Double>) {
        // Angular parts: compute the spinor spherical harmonics.
        // The four components use indices as in the note:
        //  • Component 1: Ω(θ,ϕ)_κ^(m – ½)
        //  • Component 2: Ω(θ,ϕ)_κ^(m + ½)
        //  • Component 3: Ω(θ,ϕ)_₋κ^(m – ½)
        //  • Component 4: Ω(θ,ϕ)_₋κ^(m + ½)
        let thetaPhase = Complex(cos(theta), sin(theta))
        let Y1 = sphericalHarmonics(downSpinor: false, positive: true, thetaPhase: thetaPhase, phi: phi)
        let Y2 = sphericalHarmonics(downSpinor: true, positive: true, thetaPhase: thetaPhase, phi: phi)
        let Y3 = sphericalHarmonics(downSpinor: false, positive: false, thetaPhase: thetaPhase, phi: phi)
        let Y4 = sphericalHarmonics(downSpinor: true, positive: false, thetaPhase: thetaPhase, phi: phi)
        return (Y1, Y2, Y3, Y4)
    }
    
    func waveFunction(t: Double, r: Double, theta: Double, phi: Double) -> BiSpinor {
        // Radial parts: g(ρ) and f(ρ)
        let rho = 2 * lambda * r
        let (g, f) = radial(rho: rho)
        
        // Time-dependent phase factor exp(–i ε t/ħ).
        let globalPhase = Complex<Double>.exp(self.energy * Complex(t))
        
        // Overall multiplicative factor.
        let commonFactor = globalPhase * Complex(self.normalization)  // (This factor includes normalization, time evolution and r^-1.)
        
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
        
        return BiSpinor(posiUp: posUp, posiDown: posDn, negaUp: negUp, negaDown: negDn)
    }
    
    func waveFunction(t: Double, x: Double, y: Double, z: Double) -> BiSpinor {
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
            let probability = wave.density
            if !probability.isNaN {
                totalProbability += probability * pow(HydrogenOrbital.rBohr / 100, 3)
            }
        }
    }
}
let endTime = Date()

print("Total probability in space is \(totalProbability), calculated in \(endTime.timeIntervalSince(startTime)) seconds.")

/*
let startTime = Date()
let orbital = HydrogenOrbital(z: 79, n: 24, kappa: -13, mx2: 21)
let result = orbital.waveFunction(t: 0, x: 0.4 * HydrogenOrbital.rBohr, y: 0, z: 0)
let endTime = Date()

print("Wave function is: \n\(result)\nCalculated in \(endTime.timeIntervalSince(startTime)) seconds.")
*/
