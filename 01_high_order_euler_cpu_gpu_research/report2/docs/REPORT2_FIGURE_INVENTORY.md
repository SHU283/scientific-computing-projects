# Report 2 HLLC Figure and Data Inventory

This inventory tracks figures and data for the official HLLC route. Older
HLL/Rusanov figures may remain in other folders as references, but they should
not be used as final Report 2 HLLC evidence.

## Canonical HLLC Data

### Local Smoke

Data:

```text
results/report2/local_smoke/
```

Purpose:

- quick local sanity check after code changes;
- finite/positive-state checks on small smooth and discontinuous cases;
- serial/OpenMP comparison for selected cases.

### HLLC CUDA Smoke

Data:

```text
results/hllc/cuda_checks/
```

Purpose:

- confirm the Lovelace A30 CUDA build runs;
- check WENO-2, WENO-3, and ADER-3 CUDA paths on representative smooth and discontinuous cases;
- retain the post-WENO-2-fix smoke outputs.

### Smooth 2D Convergence

Data and figure:

```text
results/hllc/validation_2d_lovelace/
report2_high_order_convergence_2d.csv
report2_high_order_convergence_2d.l1_rho.png
```

Purpose:

- confirm the HLLC route preserves the expected second/third-order smooth convergence behaviour;
- use the `200 x 200` to `400 x 400` interval as the clearest observed-order summary.

### ADER-3 CFL Diagnosis

Data:

```text
results/hllc/ader3_cfl_probe/
results/hllc/ader3_cfl012_confirm/
results/hllc/ader3_quadrant_cfl_probe/
results/hllc/ader3_discontinuous_cfl_probe/
```

Purpose:

- explain why discontinuous ADER-3 uses a more conservative CFL than WENO-2/WENO-3;
- document that CFL sensitivity is case-dependent and non-monotone, because changing CFL changes the accepted time-step sequence and can move the calculation onto different HLLC/fallback branches;
- record the ADER-3 CFL diagnosis: `0.10` fixed the checked quadrant-Riemann field but left a shock-bubble difference in the refreshed backend-consistency check;
- use the `0.13/0.14` two-case probe result, where `0.14` gave the best current compromise, before producing final figures.

## Current HLLC Evidence To Use

### Backend Consistency

Current diagnostic folder:

```text
results/hllc/backend_consistency_2d_lovelace/
```

Status:

- refreshed after the ADER-3 discontinuous CFL update;
- WENO-2/WENO-3 checks are at roundoff level;
- ADER-3 smooth Gaussian and quadrant-Riemann checks are at roundoff level;
- ADER-3 shock-bubble retains a small remaining serial/CUDA difference of about
  `2e-7`, with positive density and pressure.

Expected final evidence:

- 2D Gaussian, quadrant Riemann, and shock-bubble cases;
- WENO-2, WENO-3, and ADER-3;
- serial, 16-thread OpenMP, and CUDA;
- `100 x 100` representative fields;
- max absolute differences versus serial and positivity checks.

## Completed Performance Evidence

### Smooth Performance

Canonical data:

```text
results/hllc/performance/
```

Resolution sweep:

```text
25, 50, 100, 200, 600
```

Purpose:

- runtime scaling;
- speedup relative to serial;
- accuracy-runtime tradeoff on the smooth 2D Gaussian problem.

Timing rule:

- field CSV output is disabled inside timed runs.

### Discontinuous Performance

Canonical data:

```text
results/hllc/performance/discontinuous/
```

Resolution matrix:

```text
50, 100, 200
```

Large-grid endpoint:

```text
600 x 600
```

Purpose:

- compare backend speedup on the quadrant Riemann and shock-bubble cases;
- check whether the heavier ADER-3 discontinuous workload makes CUDA more attractive;
- keep the `600 x 600` endpoint separate because serial ADER-3 may be a long run.

### Additional Performance Studies

Canonical supporting studies under `results/hllc/performance/` include:

- OpenMP thread scaling;
- solver phase profiling;
- the accepted CUDA-only limit sweep through `5000 x 5000`;
- strict matched HLL--HLLC comparison outputs under
  `analysis/hll_vs_hllc/`.

## Figure Priority For Supervisor Updates

Use this order when a short update needs representative plots:

1. HLLC 2D convergence plot.
2. ADER-3 CFL diagnostic table/summary.
3. Canonical backend-consistency max-difference plot.
4. High-resolution 1D/2D field comparisons.
5. Smooth and discontinuous performance plots.

## Superseded Outputs

Pre-policy validation, backend metadata, and one-dimensional high-resolution
outputs are archived under
`_archive_completed_work_2026-07-19/superseded_cfl_results_2026-07-21/`.
They must not be used for the report or submission package. CFL probe folders
remain in the working tree only as documented selection evidence.
