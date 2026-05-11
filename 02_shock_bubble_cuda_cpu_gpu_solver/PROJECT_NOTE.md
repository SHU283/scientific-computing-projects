# Project Note

This is a course/research project focused on a CUDA implementation of a two-dimensional finite-volume Euler solver.

## What This Project Demonstrates

- Implemented a CPU/GPU comparison workflow for a shock-bubble benchmark.
- Used Godunov-HLL and MUSCL-Hancock finite-volume methods.
- Implemented CUDA-side field update, boundary handling, output, and runtime comparison workflow.
- Used Python scripts to post-process density and pressure fields and generate comparison figures.

## What Is Excluded

- Large generated CSV files in `data/`.
- Local build directories and temporary files.

These files can be regenerated from the source code if needed.
