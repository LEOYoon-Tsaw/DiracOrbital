//
//  model.swift
//  DiracOrbital
//
//  Created by Leo Liu on 3/23/25.
//

import Foundation
import Numerics

public struct Vector {
    let storage: InlineArray<4, Double>
    
    init(_ t: Double, _ x: Double, _ y: Double, _ z: Double) {
        storage = [t, x, y, z]
    }
    
    var abs: Double {
        sqrt(storage[0] * storage[0] - storage[1] * storage[1] - storage[2] * storage[2] - storage[3] * storage[3])
    }
}

public struct Tensor {
    let storage: InlineArray<4, InlineArray<4, Double>>
    
    init(_ storage: InlineArray<4, InlineArray<4, Double>>) {
        self.storage = storage
    }
    
    var energy: Double {
        storage[0][0]
    }
    
    var momentum: InlineArray<3, Double> {
        [storage[0][1], storage[0][2], storage[0][3]]
    }
}

public struct BispinorDerivative {
    let energy: Double
    let commonFactor: Complex<Double>
    let spherDeriv: InlineArray<4, InlineArray<3, Complex<Double>>>
    let spherCoorDeriv: InlineArray<3, InlineArray<3, Double>>
    
    func compute(waveFunction: Bispinor) async -> InlineArray<4, InlineArray<4, Complex<Double>>> {
        var matrix: InlineArray<4, InlineArray<4, Complex<Double>>> = .init(repeating: .init(repeating: 0))
        for i in 0..<4 {
            for j in 0..<3 {
                for k in 0..<3 {
                    matrix[j+1][i] += spherDeriv[i][k] * Complex(spherCoorDeriv[k][j])
                }
                matrix[j+1][i] *= -commonFactor
            }
            matrix[0][i] = Complex<Double>(0, -energy) * waveFunction.storage[i]
        }
        return matrix
    }
}

private func dot(lhs: Bispinor, rhs: Bispinor) -> Complex<Double> {
    var result = Complex<Double>(0, 0)
    for i in 0..<4 {
        result += lhs.storage[i].conjugate * rhs.storage[i]
    }
    return result
}

public struct Bispinor {
    let storage: InlineArray<4, Complex<Double>>
    
    init(_ posUp: Complex<Double>, _ posDn: Complex<Double>, _ negUp: Complex<Double>, _ negDn: Complex<Double>) {
        storage = [posUp, posDn, negUp, negDn]
    }
    init(_ storage: InlineArray<4, Complex<Double>>) {
        self.storage = storage
    }
    
    var density: Double {
        (storage[0].conjugate * storage[0] + storage[1].conjugate * storage[1] + storage[2].conjugate * storage[2] + storage[3].conjugate * storage[3]).real
    }
    
    var flow: InlineArray<3, Double> {
        [(storage[3].conjugate * storage[0] + storage[2].conjugate * storage[1] + storage[1].conjugate * storage[2] + storage[0].conjugate * storage[3]).real,
         (-storage[3].conjugate * storage[0] + storage[2].conjugate * storage[1] - storage[1].conjugate * storage[2] + storage[0].conjugate * storage[3]).imaginary,
         (storage[2].conjugate * storage[0] - storage[3].conjugate * storage[1] + storage[0].conjugate * storage[2] - storage[1].conjugate * storage[3]).real,]
    }
    
    var spin: InlineArray<3, Double> {
        [(storage[0].conjugate * storage[1] + storage[1].conjugate * storage[0] + storage[2].conjugate * storage[3] + storage[3].conjugate * storage[2]).real / 2,
         (storage[0].conjugate * storage[1] - storage[1].conjugate * storage[0] + storage[2].conjugate * storage[3] - storage[3].conjugate * storage[2]).imaginary / 2,
         (storage[0].conjugate * storage[0] - storage[1].conjugate * storage[1] + storage[2].conjugate * storage[2] - storage[3].conjugate * storage[3]).real / 2,]
    }
    
    var densityFlow: Vector {
        let flow = self.flow
        return Vector(density, flow[0], flow[1], flow[2])
    }
}

private extension Bispinor {
    static func *(lhs: Bispinor, rhs: Complex<Double>) -> Bispinor {
        return Bispinor(lhs.storage[0] * rhs, lhs.storage[1] * rhs, lhs.storage[2] * rhs, lhs.storage[3] * rhs)
    }
    
    static func *(lhs: Complex<Double>, rhs: Bispinor) -> Bispinor {
        return Bispinor(rhs.storage[0] * lhs, rhs.storage[1] * lhs, rhs.storage[2] * lhs, rhs.storage[3] * lhs)
    }
}

public extension Bispinor {
    func energyTensor(waveFDerivatives: BispinorDerivative) async -> Tensor {
        let derivative = await waveFDerivatives.compute(waveFunction: self)
        var results: InlineArray<4, InlineArray<4, Double>> = .init(repeating: .init(repeating: 0))
        for i in 0..<4 {
            let gamma2: Complex<Double> = if i % 2 == 0 { -.i } else { .i }
            let gamma3: Complex<Double> = if i % 2 == 0 { .one } else { -.one }
            let j = 3 - i
            let k = (i + 2) % 4
            results[0][0] += 4 * (storage[i].conjugate * derivative[0][i]).imaginary
            results[0][1] += (storage[i].conjugate * derivative[1][i] - derivative[1][i].conjugate * storage[i] + 2 * storage[i].conjugate * derivative[0][j]).imaginary
            results[0][2] += (storage[i].conjugate * derivative[2][i] - derivative[2][i].conjugate * storage[i] + 2 * gamma2 * storage[i].conjugate * derivative[0][j]).imaginary
            results[1][1] += (2 * storage[i].conjugate * derivative[1][j] - 2 * derivative[1][i].conjugate * storage[j]).imaginary
            results[1][2] += (storage[i].conjugate * derivative[2][j] + gamma2 * storage[i].conjugate * derivative[1][j] - derivative[2][i].conjugate * storage[j] - gamma2 * derivative[1][i].conjugate * storage[j]).imaginary
            results[1][3] += (storage[i].conjugate * derivative[3][j] + gamma3 * storage[i].conjugate * derivative[1][k] - derivative[3][i].conjugate * storage[j] - gamma3 * derivative[1][i].conjugate * storage[k]).imaginary
            results[2][2] += (2 * gamma2 * storage[i].conjugate * derivative[2][j] - 2 * gamma2 * derivative[2][i].conjugate * storage[j]).imaginary
            results[2][3] += (gamma2 * storage[i].conjugate * derivative[3][j] + gamma3 * storage[i].conjugate * derivative[2][k]).imaginary
            results[2][3] -= (gamma2 * derivative[3][i].conjugate * storage[j] - gamma3 * derivative[2][i].conjugate * storage[k]).imaginary
            results[3][3] += (2 * gamma3 * storage[i].conjugate * derivative[3][k] - 2 * gamma3 * derivative[3][i].conjugate * storage[k]).imaginary
        }
        for mu in 0..<4 {
            for nu in 0..<4 {
                if nu < mu {
                    results[mu][nu] = results[nu][mu]
                } else {
                    results[mu][nu] *= -0.25
                }
            }
        }
        return Tensor(results)
    }
}

private struct IntermediateResults {
    var zalpha: Double?
    var gamma: Double?
    var lambda: Double?
    var gfLambda: Double?
    var relativeEnergy: Double?
    var normalization: Double?
    var spinorCoefs: InlineArray<4, Complex<Double>>?
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
private func associatedLegendre(l: Int, absM: Int, cosTheta: Double) -> Double {
    if l < absM {
        return 0.0
    } else if l == absM {
        return doubleFactorial(2 * absM - 1)
    } else {
        return ((Double(2 * l - 1) * cosTheta * associatedLegendre(l: l-1, absM: absM, cosTheta: cosTheta)) - (Double(l + absM - 1) * associatedLegendre(l: l-2, absM: absM, cosTheta: cosTheta))) / Double(l - absM)
    }
}

// --- Spherical Harmonic Function ---
// Computes Y_l^m(theta, phi) = norm * P_l^m(cos(theta)) * exp(i * m * phi)
private func sphericalHarmonic(l: Int, m: Int, theta: Double, phi: Double? = nil) -> Complex<Double> {
    if l < 0 {
        return sphericalHarmonic(l: -l - 1, m: m, theta: theta, phi: phi)
    }
    let absM = abs(m)

    let sign = (m < 0 || m.isMultiple(of: 2)) ? 1.0 : -1.0
    let normFactor = sign * sqrt((Double(2 * l + 1) / (4 * Double.pi)) * factorial(l - absM) / factorial(l + absM))

    let legendreValue = associatedLegendre(l: l, absM: absM, cosTheta: cos(theta)) * Double.pow(sin(theta), absM)
    
    var result = Complex<Double>(normFactor) * Complex(legendreValue)
    if let phi {
        result *= Complex.exp(Complex(imaginary: Double(m) * phi))
    }
    
    return result
}

// --- Generalized Laguerre Polynomial ---
// Implements Lₙ^(α)(x) via a recurrence relation.
// (For n = 0: 1, for n = 1: α+1−x, and for n > 1 the recurrence is used.)
private func generalizedLaguerre(n: Int, alpha: Double, x: Double) -> Double {
    if n < 0 {
        return 0.0
    } else if n == 0 {
        return 1.0
    } else {
       return ((Double(2 * n - 1) + alpha - x) * generalizedLaguerre(n: n-1, alpha: alpha, x: x) - (Double(n - 1) + alpha) * generalizedLaguerre(n: n-2, alpha: alpha, x: x)) / Double(n)
    }
}

private func spinorCoef(a: Int, b: Int) -> Double {
    return sqrt((0.5 * Double(b + 1) + Double(a)) / Double(2 * a + 1))
}

private func sphericalHarmonicDerivatives(l: Int, m: Int, theta: Double, phi: Double) -> InlineArray<3, Complex<Double>> {
    let expFactor = Complex.exp(Complex(imaginary: Double(m) * phi))
    let y1 = sphericalHarmonic(l: l, m: m, theta: theta) * expFactor
    let y2 = sphericalHarmonic(l: l, m: m + 1, theta: theta) * expFactor
    var yDTheta = Complex(Double(m) / tan(theta)) * y1
    yDTheta += Complex<Double>.sqrt(Complex<Double>((l - m) * (l + m + 1))) * y2
    return [y1, yDTheta, .init(0, Double(m)) * y1]
}

private func rThetaPhiDXYZ(r: Double, theta: Double, phi: Double) -> InlineArray<3, InlineArray<3, Double>> {
    return [
        [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)],
        [cos(theta) * cos(phi) / r, cos(theta) * sin(phi) / r, -sin(theta) / r],
        [-sin(phi) / r / sin(theta), cos(phi) / r / sin(theta), 0.0]
    ]
}

public actor HydrogenOrbital {
    static let alpha = 0.0072973525643                     // fine structure constant
    static let c = 299792458.0                             // speed of light (m/s)
    static let hbar = 1.054571817e-34                      // reduced Planck constant (J·s)
    static let me = 9.10938356e-31                         // electron mass in kg
    static let halfPeriod = Double.pi / (mechbar * c)      // half period of global phase change
    static let mechbar = 1.0 / alpha                       // me * c / hbar
    static let rBohr = 1.0 / mechbar / alpha               // Bohr radius in m
    
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
    
    private var spinorCoefs: InlineArray<4, Complex<Double>> {
        if let spinorCoefs = intermediateResults.spinorCoefs {
            return spinorCoefs
        } else {
            let spinorCoefs: InlineArray<4, Complex<Double>> = [
                Complex(spinorCoef(a: kappa, b: -mx2)),
                Complex(-sgn(kappa) * spinorCoef(a: kappa, b: mx2)),
                Complex(0, spinorCoef(a: -kappa, b: -mx2)),
                Complex(0, sgn(kappa) * spinorCoef(a: -kappa, b: mx2))
            ]
            intermediateResults.spinorCoefs = spinorCoefs
            return spinorCoefs
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
                normalization = lambda * 0.5 * Double.pow(zalpha / gamma, 2) / Double(kappa * kappa) / Double.gamma(1 + 2 * gamma)
            }
            normalization = sqrt(normalization / commonFactor)
            intermediateResults.normalization = normalization
            return normalization
        }
    }
    
    // --- Radial functions ---
    // The radial part is split into two functions g(ρ) and f(ρ) as in the note.
    private func radial(rho: Double) -> InlineArray<2, Double> {
        let commonFactor = Double.pow(rho, gamma) * exp(-rho / 2)
        let firstTerm = rho * generalizedLaguerre(n: self.nk - 1, alpha: 2 * self.gamma + 1, x: rho)
        let lastTerm = self.gfLambda * generalizedLaguerre(n: self.nk, alpha: 2 * self.gamma - 1, x: rho)
        let gammaKappa = self.gamma - Double(self.kappa)
        let gFunc = commonFactor * (zalpha * firstTerm + gammaKappa * lastTerm)
        let fFunc = commonFactor * (zalpha * lastTerm + gammaKappa * firstTerm)
        
        return [gFunc, fFunc]
    }
    
    private func spherical(theta: Double, phi: Double) -> InlineArray<4, Complex<Double>> {
        let Y1 = spinorCoefs[0] * sphericalHarmonic(l: kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        let Y2 = spinorCoefs[1] * sphericalHarmonic(l: kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        let Y3 = spinorCoefs[2] * sphericalHarmonic(l: -kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        let Y4 = spinorCoefs[3] * sphericalHarmonic(l: -kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        return [Y1, Y2, Y3, Y4]
    }
    
    private func sphericalDerivatives(theta: Double, phi: Double) -> InlineArray<4, InlineArray<3, Complex<Double>>> {
        var o1 = sphericalHarmonicDerivatives(l: kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        var o2 = sphericalHarmonicDerivatives(l: kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        var o3 = sphericalHarmonicDerivatives(l: -kappa, m: (mx2 - 1) / 2, theta: theta, phi: phi)
        var o4 = sphericalHarmonicDerivatives(l: -kappa, m: (mx2 + 1) / 2, theta: theta, phi: phi)
        for i in 0..<3 {
            o1[i] *= spinorCoefs[0]
            o2[i] *= spinorCoefs[1]
            o3[i] *= spinorCoefs[2]
            o4[i] *= spinorCoefs[3]
        }
        return [o1, o2, o3, o4]
    }
    
    func waveFunction(t: Double, r: Double, theta: Double, phi: Double) -> Bispinor {
        // Radial parts: g(ρ) and f(ρ)
        let rho = 2 * lambda * r
        let gf = radial(rho: rho)
        
        // Time-dependent phase factor exp(–i ε t/ħ).
        let energy = Self.mechbar * Self.c * relativeEnergy
        let globalPhase = Complex<Double>.exp(Complex(imaginary: -energy * t))
        
        // Overall multiplicative factor.
        let commonFactor = globalPhase * Complex(normalization / r)  // (This factor includes normalization, time evolution and r^-1.)
        
        // Spherical parts: Y_κ^m
        var sCoefs = spherical(theta: theta, phi: phi)
        // Assemble the four components of the Dirac spinor:
        sCoefs[0] *= commonFactor * Complex(gf[0])
        sCoefs[1] *= commonFactor * Complex(gf[0])
        sCoefs[2] *= commonFactor * Complex(gf[1])
        sCoefs[3] *= commonFactor * Complex(gf[1])
        
        return Bispinor(sCoefs)
    }
    
    func waveFunction(t: Double, x: Double, y: Double, z: Double) -> Bispinor {
        let rxy = Double.hypot(x, y)
        let r = Double.hypot(rxy, z)
        let theta = Double.atan2(y: rxy, x: z)
        let phi = Double.atan2(y: y, x: x)
        return waveFunction(t: t, r: r, theta: theta, phi: phi)
    }
    
    // For derivatives
    private func gfTilde(rho: Double) -> InlineArray<4, Double> {
        let commonFactor = Double.pow(rho, gamma) * exp(-rho / 2)
        let l11 = generalizedLaguerre(n: self.nk - 1, alpha: 2 * self.gamma + 1, x: rho)
        let l12 = generalizedLaguerre(n: self.nk - 2, alpha: 2 * self.gamma + 2, x: rho)
        let l21 = generalizedLaguerre(n: self.nk, alpha: 2 * self.gamma - 1, x: rho)
        let l22 = generalizedLaguerre(n: self.nk - 1, alpha: 2 * self.gamma, x: rho)
        let gammaKappa = self.gamma - Double(self.kappa)
        let g = self.zalpha * rho * l11 + gammaKappa * self.gfLambda * l21
        let f = gammaKappa * rho * l11 + self.zalpha * self.gfLambda * l21
        let ghat = self.zalpha * (rho * l12 - l11) + gammaKappa * self.gfLambda * l22
        let fhat = gammaKappa * (rho * l12 - l11) + self.zalpha * self.gfLambda * l22
        return [
            2 * lambda * commonFactor * (((gamma - 1) / rho - 0.5) * g - ghat),
            2 * lambda * commonFactor * (((gamma - 1) / rho - 0.5) * f - fhat),
            commonFactor * g,
            commonFactor * f,
        ]
    }
    
    private func sphericalCoorDerivatives(r: Double, theta: Double, phi: Double) -> InlineArray<4, InlineArray<3, Complex<Double>>> {
        let rho = 2 * lambda * r
        let gf = gfTilde(rho: rho)
        let sphericalY = sphericalDerivatives(theta: theta, phi: phi)
        let o0: InlineArray<3, Complex<Double>> = [Complex(gf[0]) * sphericalY[0][0], Complex(gf[2]) * sphericalY[0][1], Complex(gf[2]) * sphericalY[0][2]]
        let o1: InlineArray<3, Complex<Double>> = [Complex(gf[0]) * sphericalY[1][0], Complex(gf[2]) * sphericalY[1][1], Complex(gf[2]) * sphericalY[1][2]]
        let o2: InlineArray<3, Complex<Double>> = [Complex(gf[1]) * sphericalY[2][0], Complex(gf[3]) * sphericalY[2][1], Complex(gf[3]) * sphericalY[2][2]]
        let o3: InlineArray<3, Complex<Double>> = [Complex(gf[1]) * sphericalY[3][0], Complex(gf[3]) * sphericalY[3][1], Complex(gf[3]) * sphericalY[3][2]]
        return [o0, o1, o2, o3]
    }
    
    func waveFunctionDerivatives(t: Double, r: Double, theta: Double, phi: Double) -> BispinorDerivative {
        let sphericalCoor = sphericalCoorDerivatives(r: r, theta: theta, phi: phi)
        let coorDerivatives = rThetaPhiDXYZ(r: r, theta: theta, phi: phi)
        let energy = Self.mechbar * Self.c * relativeEnergy
        let globalPhase = Complex<Double>.exp(Complex(imaginary: -energy * t))
        let commonFactor = globalPhase * Complex(normalization / r)
        
        return BispinorDerivative(energy: Self.mechbar * relativeEnergy, commonFactor: commonFactor, spherDeriv: sphericalCoor, spherCoorDeriv: coorDerivatives)
    }
    
    func waveFunctionDerivatives(t: Double, x: Double, y: Double, z: Double) -> BispinorDerivative {
        let rxy = Double.hypot(x, y)
        let r = Double.hypot(rxy, z)
        let theta = Double.atan2(y: rxy, x: z)
        let phi = Double.atan2(y: y, x: x)
        return waveFunctionDerivatives(t: t, r: r, theta: theta, phi: phi)
    }
}
