# Project Note

This is a course/research project focused on OpenMP and MPI parallelization for a two-dimensional compressible Euler solver.

## What This Project Demonstrates

- Extended a serial finite-volume solver to OpenMP and MPI versions.
- Implemented domain decomposition, halo exchange, and rank-0 output gathering for the MPI version.
- Used MUSCL-Hancock reconstruction, HLL flux, Strang splitting, and CFL-based time stepping.
- Performed runtime, speedup, and parallel-efficiency analysis with Python post-processing scripts.
