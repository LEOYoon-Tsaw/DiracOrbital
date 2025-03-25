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
        
        var numberOfTasks = 16
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
        var tasks = [Task<[(Double, Double)], Never>]()
        
        for i in 0..<numberOfTasks {
            let task = Task {
                let orbital = HydrogenOrbital(z: z, n: n, kappa: kappa, mx2: mx2)
                var results = [(Double, Double)]()
                
                for r in 0...totalRRuns {
                    let rAct = Double(r) * deltaR * HydrogenOrbital.rBohr
                    for theta in (i * totalThetaRuns / numberOfTasks)..<((i+1) * totalThetaRuns / numberOfTasks) {
                        let thetaAct = Double(theta) * deltaTheta * Double.pi
                        let wave = await orbital.waveFunction(t: 0, r: Double(r) * deltaR * HydrogenOrbital.rBohr, theta: Double(theta) * deltaTheta * Double.pi, phi: 0)
                        var probability = 0.0
                        for component in wave {
                            probability += (component.conjugate * component).real
                        }
                        if !probability.isNaN {
                            let density = probability
                            probability *= 2 * Double.pi * rAct * rAct * sin(thetaAct) * dr * dTheta
                            results.append((probability, density))
                        }
                    }
                }
                
                return results
            }
            
            tasks.append(task)
        }
        
        var allResults: [(Double, Double)] = []
        
        for task in tasks {
            let results = await task.value
            allResults.append(contentsOf: results)
        }
        
        allResults.sort { $0.1 > $1.1 }
        let benchmarks = [0.25, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.97, 0.98, 0.99]
        var sum = 0.0
        var i = 0
        var percentiles = [Double: Double]()
        for (prob, density) in allResults {
            sum += prob
            while i < benchmarks.count && sum >= benchmarks[i] {
                percentiles[sum] = density
                i += 1
            }
        }
        
        let endTime = Date()
        
        for (percentile, density) in percentiles.sorted(by: { $0.key < $1.key }) {
            print(String(format: "%.1f%%: density: %.5e", percentile * 100, density))
        }
        print("Covers \(String(format: "%.3f%%", sum * 100)) of entire space. Calculated in \(endTime.timeIntervalSince(startTime)) seconds.")
    }
}
