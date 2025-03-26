## Example of Usage
In terminal: `(base) > DiracOrbital -Z 79 -n 4 -k -1 -m -1`

Output:
```
Starting task:
Z=79, n=4, k=-1 (l=1, j=1/2), m=-1/2
r ranging from 0 to 1.0 Bohr radius, on 5000 × 1000 (r×θ) grid, running on 10 threads.
25.0%: density: 3.75474e+31, at 0.274 Bohr radius, 0.012π theta
50.0%: density: 2.98915e+31, at 0.305 Bohr radius, 0.008π theta
55.0%: density: 2.75473e+31, at 0.312 Bohr radius, 0.000π theta
60.0%: density: 2.49902e+31, at 0.319 Bohr radius, 0.001π theta
65.0%: density: 2.22127e+31, at 0.328 Bohr radius, 0.003π theta
70.0%: density: 1.92479e+31, at 0.338 Bohr radius, 0.000π theta
75.0%: density: 1.61425e+31, at 0.348 Bohr radius, 0.003π theta
80.0%: density: 1.29086e+31, at 0.361 Bohr radius, 0.002π theta
85.0%: density: 9.58670e+30, at 0.377 Bohr radius, 0.031π theta
90.0%: density: 6.25823e+30, at 0.398 Bohr radius, 0.001π theta
95.0%: density: 2.96494e+30, at 0.431 Bohr radius, 0.000π theta
97.0%: density: 1.70245e+30, at 0.454 Bohr radius, 0.018π theta
98.0%: density: 1.09025e+30, at 0.472 Bohr radius, 0.001π theta
99.0%: density: 5.11989e+29, at 0.501 Bohr radius, 0.003π theta
Covers 100.000% of entire space. Calculated in 6.264032006263733 seconds.
```

## All Parameters

In terminal: `(base) > DiracOrbital -h`

Output:
```
Possible arguments: -Z, -n, -k (kappa, Dirac quantum number), -m (2 × m, integer), and --deltaR (float in Bohr radius), --deltaTheta (float, 1 = π), --totalR (float in Bohr radius), --numberOfTasks
```
