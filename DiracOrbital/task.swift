//
//  task.swift
//  DiracOrbital
//
//  Created by Leo Liu on 3/24/25.
//

import Foundation

enum Mode: String {
    case tasks = "--numberOfTasks", z = "-Z", n = "-n", kappa = "-k", m = "-m"
    case deltaR = "--deltaR", deltaTheta = "--deltaTheta", totalR = "--totalR"
    case help = "--help", h = "-h"
}

@main
struct Main {
    static func main() async {
        
        var numberOfTasks = 10
        var z = 79
        var n = 4
        var kappa = -1
        var mx2 = 1
        var deltaR = 0.0002
        var deltaTheta = 0.001
        var totalR = 1.0
        
        var mode: Mode?
        let args = CommandLine.arguments
        for arg in args[1..<args.count] {
            if let _mode = mode {
                switch _mode {
                case .tasks:
                    numberOfTasks = Int(arg) ?? numberOfTasks
                case .z:
                    z = Int(arg) ?? z
                case .n:
                    n = Int(arg) ?? n
                case .kappa:
                    kappa = Int(arg) ?? kappa
                case .m:
                    mx2 = Int(arg) ?? mx2
                case .deltaR:
                    deltaR = Double(arg) ?? deltaR
                case .deltaTheta:
                    deltaTheta = Double(arg) ?? deltaTheta
                case .totalR:
                    totalR = Double(arg) ?? totalR
                case .help, .h:
                    break
                }
                mode = nil
            } else {
                mode = Mode(rawValue: arg)
                if mode == .help || mode == .h {
                    print("Possible arguments: -Z, -n, -k (kappa, Dirac quantum number), -m (2 × m, integer), and --deltaR (float in Bohr radius), --deltaTheta (float, 1 = π), --totalR (float in Bohr radius), --numberOfTasks")
                    return
                }
            }
        }
        
        let dr = deltaR * HydrogenOrbital.rBohr
        let dTheta = deltaTheta * Double.pi
        let totalThetaRuns = Int(1 / deltaTheta) + 1
        let totalRRuns = Int(totalR / deltaR)
        
        let jx2 = 2 * abs(kappa) - 1
        let l = kappa > 0 ? kappa + 1 : -kappa
        print("Starting task:\nZ=\(z), n=\(n), k=\(kappa) (l=\(l), j=\(jx2)/2), m=\(mx2)/2\nr ranging from 0 to \(totalR) Bohr radius, on \(totalRRuns) × \(totalThetaRuns-1) (r×θ) grid, running on \(numberOfTasks) threads.")
        
        let startTime = Date()
        var tasks = [Task<[(Double, Double, Double, Double, Double, Double, Double)], Never>]()
        
        for i in 0..<numberOfTasks {
            let task = Task {
                let orbital = HydrogenOrbital(z: z, n: n, kappa: kappa, mx2: mx2)
                var results = [(Double, Double, Double, Double, Double, Double, Double)]()
                
                for r in 0...totalRRuns {
                    let rAct = Double(r) * deltaR * HydrogenOrbital.rBohr
                    for theta in (i * totalThetaRuns / numberOfTasks)..<((i+1) * totalThetaRuns / numberOfTasks) {
                        let thetaAct = Double(theta) * deltaTheta * Double.pi
                        let x = rAct * sin(thetaAct)
                        let dThetaDist = rAct * dTheta
                        let phiDist = 2 * Double.pi * x
                        let wave = await orbital.waveFunction(t: 0, r: rAct, theta: thetaAct, phi: 0)
                        let waveDeriv = await orbital.waveFunctionDerivatives(t: 0, r: rAct, theta: thetaAct, phi: 0)
                        let orbitalMom = x * dot(lhs: wave, rhs: waveDeriv[1]).imaginary * phiDist * dThetaDist * dr
                        let density = wave.density
                        let flow = wave.flow
                        let spin = wave.spin
                        if !density.isNaN && !density.isZero {
                            let probability = density * phiDist * dThetaDist * dr
                            let magnetic = flow[1] * x * phiDist * dThetaDist * dr
                            let actSpin = spin[2] * phiDist * dThetaDist * dr
                            results.append((probability, density, rAct / HydrogenOrbital.rBohr, thetaAct, magnetic * HydrogenOrbital.mechbar, actSpin, orbitalMom))
                        }
                    }
                }
                
                return results
            }
            
            tasks.append(task)
        }
        
        var allResults: [(Double, Double, Double, Double, Double, Double, Double)] = []
        
        for task in tasks {
            let results = await task.value
            allResults.append(contentsOf: results)
        }
        
        allResults.sort { $0.1 > $1.1 }
        let benchmarks = [0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.97, 0.98, 0.99]
        var sum = 0.0
        var magnSum = 0.0
        var spinSum = 0.0
        var orbitalSum = 0.0
        var i = 0
        var maxR = 0.0
        var maxTheta = 0.0
        var percentiles = [Double: (Double, Double, Double, Double, Double, Double)]()
        for (prob, density, r, theta, magnetic, spin, orbital) in allResults {
            sum += prob
            magnSum += magnetic
            spinSum += spin
            if !orbital.isNaN {
                orbitalSum += orbital
            }
            if r > maxR {
                maxR = r
                maxTheta = theta
            }
            while i < benchmarks.count && sum >= benchmarks[i] {
                percentiles[sum] = (density, maxR, maxTheta, magnSum, spinSum, orbitalSum)
                i += 1
            }
        }
        
        let endTime = Date()
        
        for (percentile, (density, r, theta, magnetic, spin, orbital)) in percentiles.sorted(by: { $0.key < $1.key }) {
            print(String(format: "%.1f%%: density: %.5e, at %.3f Bohr radius, %.3fπ theta; covered magnetic momentum: %.3fμB, spin: %.3fħ, orbital: %.3fħ", percentile * 100, density, r, theta < 0 ? 1 + theta / Double.pi : theta / Double.pi, magnetic, spin, orbital))
        }
        print("Covers \(String(format: "%.3f%%", sum * 100)) of entire space, cumulating \(String(format: "%.3f", magnSum))μB magnetic momentum, \(String(format: "%.3f", spinSum))ħ spin, \(String(format: "%.3f", orbitalSum))ħ orbital. Calculated in \(String(format: "%.2f", endTime.timeIntervalSince(startTime))) seconds.")
    }
}
