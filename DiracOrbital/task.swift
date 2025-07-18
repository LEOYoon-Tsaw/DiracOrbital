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
    case angularMomentum = "--angularMomentum"
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
        var angularMomentum: Bool = false
        
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
                default:
                    break
                }
                mode = nil
            } else {
                mode = Mode(rawValue: arg)
                switch mode {
                case .angularMomentum:
                    angularMomentum = true
                case .help, .h:
                    print("Possible arguments: -Z (int), -n (int), -k (kappa, Dirac quantum number), -m (2 × m, integer), and --deltaR (float in Bohr radius), --deltaTheta (float, 1 = π), --totalR (float in Bohr radius), --numberOfTasks (int), and --angularMomentum")
                    return
                default:
                    break
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
        var tasks = [Task<[InlineArray<6, Double>], Never>]()
        
        for i in 0..<numberOfTasks {
            let task = Task {
                let orbital = HydrogenOrbital(z: z, n: n, kappa: kappa, mx2: mx2)
                var results = [InlineArray<6, Double>]()
                
                for r in 0...totalRRuns {
                    let rAct = Double(r) * deltaR * HydrogenOrbital.rBohr
                    for theta in (i * totalThetaRuns / numberOfTasks)..<((i+1) * totalThetaRuns / numberOfTasks) {
                        let thetaAct = Double(theta) * deltaTheta * Double.pi
                        let x = rAct * sin(thetaAct)
                        let dThetaDist = rAct * dTheta
                        let phiDist = 2 * Double.pi * x
                        let wave = await orbital.waveFunction(t: 0, r: rAct, theta: thetaAct, phi: 0)
                        let density = wave.density
                        let flow = wave.flow
                        if !density.isNaN && !density.isZero {
                            let probability = density * phiDist * dThetaDist * dr
                            let magnetic = -flow[1] * x * phiDist * dThetaDist * dr
                            var angular: Double = .nan
                            if angularMomentum {
                                let waveDeriv = await orbital.waveFunctionDerivatives(t: 0, r: rAct, theta: thetaAct, phi: 0)
                                let energeStressTensor = await wave.energyTensor(waveFDerivatives: waveDeriv)
                                angular = energeStressTensor.momentum[1] * x * phiDist * dThetaDist * dr
                            }
                            results.append([probability, density, rAct / HydrogenOrbital.rBohr, thetaAct, magnetic * HydrogenOrbital.mechbar, angular])
                        }
                    }
                }
                
                return results
            }
            
            tasks.append(task)
        }
        
        var allResults: [InlineArray<6, Double>] = []
        
        for task in tasks {
            let results = await task.value
            allResults.append(contentsOf: results)
        }
        
        allResults.sort { $0[1] > $1[1] } // index 1 is prob density
        let benchmarks = [0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.97, 0.98, 0.99]
        var sum = 0.0
        var magnSum = 0.0
        var angSum = 0.0
        var i = 0
        var maxR = 0.0
        var maxTheta = 0.0
        var percentiles = [Double: InlineArray<5, Double>]()
        for values in allResults {
            sum += values[0] // prob
            magnSum += values[4] // magnetic
            if !values[5].isNaN { // angular
                angSum += values[5]
            }
            if values[2] > maxR { // r
                maxR = values[2]
                maxTheta = values[3] // theta
            }
            while i < benchmarks.count && sum >= benchmarks[i] {
                percentiles[sum] = [values[0], maxR, maxTheta, magnSum, angSum]
                i += 1
            }
        }
        
        let endTime = Date()
        
        for (percentile, values) in percentiles.sorted(by: { $0.key < $1.key }) {
            var line = String(format: "%.1f%% coverage -> density: %.5e, at %.3f rB, θ = %.3fπ; covers magnetic moment: %.3fμB", percentile * 100, values[0], values[1], values[2] < 0 ? 1 + values[2] / Double.pi : values[2] / Double.pi, values[3])
            if angularMomentum {
                line += String(format: ", angular momentum: %.3fħ", values[4])
            }
            line += "."
            print(line)
        }
        var line = "Covers \(String(format: "%.3f%%", sum * 100)) probability, accumulating \(String(format: "%.3f", magnSum))μB magnetic moment"
        if angularMomentum {
            line += ", \(String(format: "%.3f", angSum))ħ angular momemtum"
        }
        line += ". Calculated in \(String(format: "%.2f", endTime.timeIntervalSince(startTime))) seconds."
        print(line)
    }
}
