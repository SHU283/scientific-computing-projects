# ADER-3 CFL-selection evidence

This directory contains the probe data used to select the ADER-3
discontinuous CFL value.

## Final policy

- Smooth ADER-3: CFL 0.25.
- Discontinuous ADER-3: CFL 0.14.
- Every serial/OpenMP/CUDA comparison at a formal data point uses the same CFL,
  final time, boundary condition, and grid for that row.

## Included probe groups

- `ader3_cfl012_confirm/`: focused confirmation at CFL 0.12.
- `ader3_cfl_probe/`: shock-bubble probe records around the sensitive range.
- `ader3_discontinuous_cfl_probe/`: paired quadrant-Riemann and shock-bubble
  summaries, including the CFL 0.13 and 0.14 comparison.
- `ader3_quadrant_cfl_probe/`: broader quadrant-Riemann probe fields and logs.

The summary CSVs, logs, and exported fields show backend differences and
admissibility at each tested CFL. The main performance and backend-consistency
experiments use the selected CFL and are stored under `results/hllc/`.
