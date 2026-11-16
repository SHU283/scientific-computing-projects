# HLLC phase-profile audit

- Raw rows: 36
- Expected rows: 36
- Audit findings: 0

| n | method | backend | repeats | median total (s) | range / median | accounted | process overhead (s) |
|---:|:---|:---|---:|---:|---:|---:|---:|
| 100 | WENO2 | serial | 3 | 0.702894 | 0.94% | 98.49% | 0.01416 |
| 100 | WENO2 | openmp | 3 | 0.100226 | 4.70% | 95.24% | 0.01098 |
| 100 | WENO2 | cuda | 3 | 0.0159133 | 7.56% | 97.66% | 0.5404 |
| 100 | WENO3 | serial | 3 | 2.38888 | 5.10% | 98.51% | 0.02139 |
| 100 | WENO3 | openmp | 3 | 4.76007 | 40.90% | 99.32% | 0.01783 |
| 100 | WENO3 | cuda | 3 | 0.0379426 | 3.88% | 99.19% | 0.5557 |
| 100 | ADER3 | serial | 3 | 50.2444 | 4.33% | 100.00% | 0.02213 |
| 100 | ADER3 | openmp | 3 | 8.74594 | 1.03% | 99.99% | 0.02895 |
| 100 | ADER3 | cuda | 3 | 0.0747606 | 112.72% | 99.56% | 0.6048 |
| 600 | WENO2 | serial | 1 | 159.788 | 0.00% | 98.55% | 0.3438 |
| 600 | WENO2 | openmp | 1 | 19.5098 | 0.00% | 94.38% | 0.1824 |
| 600 | WENO2 | cuda | 1 | 0.862241 | 0.00% | 99.04% | 0.7249 |
| 600 | WENO3 | serial | 1 | 424.145 | 0.00% | 98.94% | 0.2722 |
| 600 | WENO3 | openmp | 1 | 61.726 | 0.00% | 92.59% | 0.2733 |
| 600 | WENO3 | cuda | 1 | 2.24114 | 0.00% | 99.65% | 0.7437 |
| 600 | ADER3 | serial | 1 | 10171.6 | 0.00% | 100.00% | 0.2982 |
| 600 | ADER3 | openmp | 1 | 907.777 | 0.00% | 100.00% | 0.314 |
| 600 | ADER3 | cuda | 1 | 10.0399 | 0.00% | 99.93% | 0.7281 |

## Findings

- No completeness or timing-consistency problems detected.
