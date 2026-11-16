# ADER-3 CFL-selection evidence

This directory preserves the selected probe data used to justify the final
ADER-3 discontinuous CFL value reported in the dissertation. It is supporting
experimental history, not a second canonical result route.

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

The files include summary CSVs, console logs, and exported fields where they
are needed to inspect backend differences and admissibility. They document why
the practical compromise CFL 0.14 was selected. They must not be mixed into
the final performance or backend-consistency tables, whose canonical files are
under `results/hllc/`.
