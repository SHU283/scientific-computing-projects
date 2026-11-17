# Report 1 Code Submission Package

This folder contains the source code, scripts, settings files, and selected generated outputs used to support the results presented in Project Report 1.

## Contents

- `CMakeLists.txt`, `include/`, `src/`: C++/CUDA source code and headers for the finite-volume Euler solver.
- `scripts/`: Python scripts used to run experiments, compare outputs, summarize timings, and generate figures.
- `results/stageA/`: validation, convergence, and representative field outputs.
- `results/stageB/`: serial/OpenMP/CUDA backend-consistency outputs.
- `results/stageC/`: runtime, speedup, OpenMP thread-scaling, and accuracy--runtime outputs.
- `results/stageD/`: direct method-comparison figures, difference maps, and slice plots.

The final timing results used in the written report are the Lovelace A30 files in `results/stageC/`. The all-backend 2D Gaussian runtime-scaling, speedup, and accuracy--runtime files use the `lovelace_a30_noio` naming convention, because CSV field-output time is excluded from those timings. The Lovelace OpenMP thread-scaling and discontinuous runtime checks are also included. The all-backend OpenMP timings use a fixed configuration (`OMP_NUM_THREADS=16`, `OMP_DYNAMIC=FALSE`) so that they match the written report; the separate thread-scaling study intentionally varies `OMP_NUM_THREADS`.

## Build

From a Linux or WSL shell:

```bash
cmake -S . -B build
cmake --build build -j
```

CUDA build, if a CUDA compiler and GPU environment are available:

```bash
cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON
cmake --build build-cuda -j
```

## Example Runs

```bash
./build/final_project --study=convergence --problem=gaussian_advection_2d --method=muscl
./build/final_project --backend=openmp --problem=gaussian_advection_2d --method=weno2 --nx=100 --ny=100
./build-cuda/final_project --backend=cuda --problem=gaussian_advection_2d --method=weno2 --nx=100 --ny=100
```

The Python scripts in `scripts/` provide the report-specific batch runs and plotting commands.

## Methods

The main Report 1 results use Godunov, MUSCL-Hancock, and WENO-2.
Preliminary WENO-3 and ADER-3 implementations continue in Report 2.
