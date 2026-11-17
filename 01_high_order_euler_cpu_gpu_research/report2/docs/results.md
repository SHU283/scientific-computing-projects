# Report 2 numerical and performance results

The experiments use the HLLC solver with serial, OpenMP, and CUDA backends.
The HLL/Rusanov comparison is in `analysis/hll_vs_hllc/`.

## Numerical accuracy and backend agreement

The common CFL settings are 0.40 for smooth WENO-2/WENO-3 and 0.25 for
smooth ADER-3; discontinuous cases use 0.20 and 0.14 respectively.
The CFL-aligned 2D convergence campaign recorded finest-interval L1 density
orders of 2.0067, 2.9981, and 2.9980 for WENO-2, WENO-3, and ADER-3.

All 27 backend-check fields were finite with positive density and pressure.
OpenMP matched the exported serial fields exactly. CUDA differences were at
roundoff or near-roundoff level except for the branch-sensitive shock-bubble
ADER-3 case, whose maximum energy difference was approximately `2.02e-7`.
The [HLLC notes](REPORT2_HLLC_ROUTE.md) explain the CFL sensitivity and the
choice of 0.14 for discontinuous ADER-3.

The convergence and field datasets are under `results/hllc/validation_1d/`,
`validation_2d_lovelace/`, and `cfl_aligned_validation/`; the extension
through 1600 is under `results/hllc/extended_convergence_1600/`.
Backend differences are in `results/hllc/backend_consistency_2d_lovelace/`.

## Gaussian advection at 600 × 600

The source table is
`results/hllc/performance/runtime_accuracy_2d_report2_hllc_lovelace_a30.csv`.
Field CSV output is excluded from these runtime measurements.

| Method | L1 density error | Serial (s) | OpenMP (s) | CUDA (s) | OpenMP speedup | CUDA speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| WENO-2 | 4.70263e-7 | 79.85748 | 19.56291 | 2.07785 | 4.0821 | 38.4327 |
| WENO-3 | 1.72624e-8 | 227.30518 | 37.13473 | 2.71128 | 6.1211 | 83.8369 |
| ADER-3 | 1.68102e-8 | 4118.48147 | 277.20461 | 10.50597 | 14.8572 | 392.0136 |

The WENO-3 and ADER-3 density errors are lower than WENO-2 by factors of
27.2420 and 27.9749 for this case. Discontinuous runtime, OpenMP thread
scaling, phase profiles, and CUDA scaling through 5000 × 5000 are grouped
under `results/hllc/performance/`.

## WENO-2 OpenMP correction

The July 2026 correction parallelized the CFL maximum reduction and both
SSPRK2 state-update stages. The performance tables include the corrected
WENO-2 OpenMP rows. Serial/OpenMP field comparisons at 100 × 100 gave identical
exported values for Gaussian advection, quadrant Riemann, and shock-bubble.

At 600 × 600, WENO-2 OpenMP took 95.18 s for quadrant Riemann and 173.33 s
for shock-bubble. The 16-thread scaling speedup was 4.204. In the phase profile,
the CFL cost fell from 2.99 s to 0.22 s; the spatial operator accounted for
59.5% of the corrected OpenMP solver time and ran 12.20 times faster than serial.

The [refresh results](../results/hllc/performance/weno2_openmp_refresh/VALIDATION.md)
describe timing variability. [Phase-profiling notes](phase_profiling_2d.md)
define the measurement categories and diagnostic synchronizations.

## Recorded regression checks

The package validation on 2026-08-05 recorded six successful builds:
HLLC and HLL/Rusanov, each in serial-only, OpenMP, and CUDA configurations.
A 68-case smoke matrix covered WENO-2, WENO-3, ADER-3, smooth and
discontinuous 2D problems, and representative 1D cases.

The backend-alignment run covered all nine HLLC problem-method combinations
at 100 × 100. Serial and 16-thread OpenMP used identical step counts and
matched all exported decimal field values (`max_delta = 0`).

Run the local regression suite from `code/hllc_solver/` after building:

```bash
python3 scripts/run_report2_local_smoke.py --build-dir build-report2
```

The backend-field comparison, CFL-policy checks, and phase-profile analysis
scripts remain available in the same `scripts/` directory.
