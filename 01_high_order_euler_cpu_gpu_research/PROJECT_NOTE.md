# Project Note

This is an ongoing MPhil Scientific Computing project on high-order finite-volume methods for the compressible Euler equations.

## What This Project Demonstrates

- Built a unified experimental framework for comparing finite-volume methods.
- Implemented and tested Godunov-HLL, MUSCL-Hancock, and WENO-2 style methods.
- Integrated serial CPU, OpenMP, and CUDA backend paths for runtime and consistency analysis.
- Used Python scripts for convergence studies, field visualization, and performance summaries.

## Current Status

The project is still part of an ongoing thesis workflow. Godunov-HLL, MUSCL-Hancock, WENO-2, serial CPU, OpenMP, and CUDA-related components are the main completed or actively tested parts. Some planned higher-order extensions may be experimental.

## What Is Excluded

- Large generated CSV/result files in `results/` and `data/`.
- Local build directories.
- Historical writing notes and archived drafts.

These files are excluded to keep the repository focused on source code, representative outputs, and reproducible scripts.
