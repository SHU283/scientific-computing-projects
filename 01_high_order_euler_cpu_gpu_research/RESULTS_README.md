# Results Guide

This folder keeps a curated subset of outputs from the high-order Euler solver project. The full experimental workspace contains many more intermediate CSV files and plots, while this repository keeps only representative results and figures.

## Included Results

- `results/stageA/`: representative accuracy tables, 1D/2D benchmark fields, and selected raw CSV outputs.
- `results/stageB/`: serial/OpenMP/CUDA consistency summaries and selected backend-comparison fields.
- `results/stageC/`: compact runtime, OpenMP thread-scaling, accuracy-runtime, and backend-speedup summaries.
- `results/stageD/`: method-comparison figures for Godunov-HLL, MUSCL-Hancock, and WENO-2.

## Why These Files Are Kept

- They cover both smooth accuracy tests and discontinuous shock-dominated tests.
- They include both visual field outputs and raw CSV data for representative cases.
- They document method comparisons, CPU/OpenMP/CUDA consistency, runtime, and speedup.
- They avoid uploading every intermediate resolution/backend file, keeping the repository readable.

## Recommended Figures to Inspect

- `results/stageA/shock_bubble_2d_*_t020.rho.png`
- `results/stageA/shock_bubble_2d_*_t020.p.png`
- `results/stageB/backend_consistency_report1_summary.max_delta.png`
- `results/stageC/runtime_scaling_2d_report1_all_backends.best_seconds.speedup.png`
- `results/stageC/openmp_thread_scaling_report1.best_seconds.png`
- `results/stageD/shock_bubble_2d_godunov_t020__vs__shock_bubble_2d_muscl_t020__vs__shock_bubble_2d_weno2_t020.rho.png`

The selected CSV files are kept so the representative plots can be regenerated or checked numerically.
