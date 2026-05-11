# Results Guide

This folder includes representative outputs from the OpenMP/MPI shock-bubble solver project.

## Included Results

- `figs/`: field plots, OpenMP/MPI scaling plots, efficiency plots, and difference plots.
- `sample_data/openmp_final_fields/`: selected final OpenMP field CSV outputs.
- `sample_data/mpi_scaling_fields/`: selected MPI field CSV outputs for several process counts.
- `sample_data/benchmark_logs/`: strong-scaling, weak-scaling, and high-resolution MPI benchmark logs.

## Recommended Figures to Inspect

- `figs/omp_np1_np4_p_final_side_by_side.png`
- `figs/mpi_speedup.png`
- `figs/mpi_efficiency.png`
- `figs/highres_mpi_speedup.png`
- `figs/weak_runtime_per_step.png`
- `figs/rho_2000.png`

The included logs and CSV files provide the numerical basis for the scaling plots.
