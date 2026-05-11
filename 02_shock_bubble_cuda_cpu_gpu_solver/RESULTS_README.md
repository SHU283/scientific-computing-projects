# Results Guide

This folder includes representative outputs from the CUDA/CPU shock-bubble solver project.

## Included Results

- `project_figs/`: rendered pressure/density fields, CPU-GPU difference plots, lineouts, and runtime comparisons.
- `sample_data/first_order_cpu/`: selected raw CSV snapshots from the first-order CPU run.
- `sample_data/first_order_gpu/`: selected raw CSV snapshots from the first-order GPU run.
- `sample_data/second_order_cpu/`: selected raw CSV snapshots from the second-order CPU run.
- `sample_data/second_order_gpu/`: selected raw CSV snapshots from the second-order GPU run.

## Recommended Figures to Inspect

- `project_figs/rho_triple_cpu1_gpu1_absdiff_step500.png`
- `project_figs/p_triple_cpu2_gpu2_absdiff_step500.png`
- `project_figs/rho_max_abs_error_vs_step_cpu1_gpu1.png`
- `project_figs/runtime_total_cpu_vs_gpu.png`
- `project_figs/timing_breakdown.png`

Only selected CSV snapshots are included to keep the repository reviewable while still preserving raw data examples.
