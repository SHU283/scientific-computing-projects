# HLLC 2D phase-profile summary

Times are medians where repeats are available. Speedups use the profiled serial run
for the same method and resolution. CUDA phase timings include diagnostic boundary
synchronizations and are therefore not replacements for the canonical no-I/O timings.

| n | method | backend | total (s) | speedup | dominant phase | share |
|---:|:---|:---|---:|---:|:---|---:|
| 100 | ADER3 | serial | 50.24 | 1.000 | spatial_operator | 99.7% |
| 100 | ADER3 | openmp | 8.746 | 5.745 | spatial_operator | 81.6% |
| 100 | ADER3 | cuda | 0.07476 | 672.071 | spatial_operator | 77.8% |
| 100 | WENO2 | serial | 0.7029 | 1.000 | spatial_operator | 88.8% |
| 100 | WENO2 | openmp | 0.1002 | 7.013 | spatial_operator | 61.2% |
| 100 | WENO2 | cuda | 0.01591 | 44.170 | spatial_operator | 62.0% |
| 100 | WENO3 | serial | 2.389 | 1.000 | spatial_operator | 90.5% |
| 100 | WENO3 | openmp | 4.76 | 0.502 | spatial_operator | 69.0% |
| 100 | WENO3 | cuda | 0.03794 | 62.960 | spatial_operator | 68.1% |
| 600 | ADER3 | serial | 1.017e+04 | 1.000 | spatial_operator | 99.7% |
| 600 | ADER3 | openmp | 907.8 | 11.205 | spatial_operator | 98.4% |
| 600 | ADER3 | cuda | 10.04 | 1013.116 | spatial_operator | 96.3% |
| 600 | WENO2 | serial | 159.8 | 1.000 | spatial_operator | 88.6% |
| 600 | WENO2 | openmp | 19.51 | 8.190 | spatial_operator | 59.5% |
| 600 | WENO2 | cuda | 0.8622 | 185.317 | spatial_operator | 80.5% |
| 600 | WENO3 | serial | 424.1 | 1.000 | spatial_operator | 91.0% |
| 600 | WENO3 | openmp | 61.73 | 6.871 | spatial_operator | 59.5% |
| 600 | WENO3 | cuda | 2.241 | 189.254 | spatial_operator | 83.8% |

## Large-grid phase composition (600 x 600)

`setup/workspace/transfer` combines initialization, temporary flux-array allocation,
and explicit transfer timing. It is useful as a cost fraction, but its CPU and CUDA
implementation boundaries are not identical, so it is excluded from phase-speedup claims.

| method | backend | setup/workspace/transfer | CFL | spatial | update/validity | other |
|:---|:---|---:|---:|---:|---:|---:|
| WENO2 | serial | 5.0% | 2.0% | 88.6% | 2.8% | 1.4% |
| WENO2 | openmp | 27.6% | 1.1% | 59.5% | 6.2% | 5.6% |
| WENO2 | cuda | 0.5% | 9.2% | 80.5% | 8.8% | 1.0% |
| WENO3 | serial | 3.8% | 1.3% | 91.0% | 2.8% | 1.1% |
| WENO3 | openmp | 27.5% | 0.6% | 59.5% | 5.0% | 7.4% |
| WENO3 | cuda | 0.2% | 6.3% | 83.8% | 9.4% | 0.3% |
| ADER3 | serial | 0.0% | 0.1% | 99.7% | 0.2% | 0.0% |
| ADER3 | openmp | 0.4% | 0.5% | 98.4% | 0.7% | 0.0% |
| ADER3 | cuda | 0.0% | 2.6% | 96.3% | 1.0% | 0.1% |

## Comparable large-grid phase speedups

Ratios below use the profiled serial phase for the same method. They diagnose where
acceleration is obtained; they are not substitutes for canonical no-I/O speedups.

| method | backend | CFL | spatial operator | update/validity |
|:---|:---|---:|---:|---:|
| WENO2 | openmp | 14.72x | 12.20x | 3.75x |
| WENO2 | cuda | 40.84x | 204.06x | 59.48x |
| WENO3 | openmp | 14.45x | 10.52x | 3.83x |
| WENO3 | cuda | 39.72x | 205.64x | 55.91x |
| ADER3 | openmp | 2.05x | 11.36x | 2.41x |
| ADER3 | cuda | 37.19x | 1049.00x | 158.52x |

## Interpretation and data-quality notes

- The spatial operator is the dominant phase for every 600 x 600 run. Its share remains especially high for ADER-3, so further ADER-3 optimization should start with predictor/flux kernels rather than peripheral work.
- OpenMP reduces the spatial cost strongly, but WENO-2 and WENO-3 retain substantial workspace, CFL, and update overhead. This explains why their end-to-end OpenMP gains are smaller than their spatial-operator gains.
- CUDA also leaves the spatial operator dominant at large resolution. Kernel and memory-access optimization therefore remain more valuable than optimizing final diagnostics.
- The 600 x 600 rows have one repeat each, so their phase fractions are diagnostic rather than uncertainty estimates.
- Small-grid groups with at least 20% range/median variation are:
  - ADER3 cuda at 100 x 100: 3 repeats, 112.7% range/median.
  - WENO3 openmp at 100 x 100: 3 repeats, 40.9% range/median.
  These cases should illustrate startup/profiling sensitivity, not supply final performance claims.
