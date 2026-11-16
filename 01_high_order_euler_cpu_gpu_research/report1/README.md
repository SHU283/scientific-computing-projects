# Report 1 Code Package

This is a curated code package for the Project Report 1 submission.  It contains
the finite-volume solver, experiment scripts, selected outputs, and report source
needed to support the results presented in the written report.

Main implementation focus:

1. Godunov, MUSCL-Hancock, and WENO-2 finite-volume solvers
2. One- and two-dimensional Euler/advection test problems
3. Serial, OpenMP, and CUDA backend checks
4. Validation, convergence, runtime, and plotting scripts
5. Preliminary WENO-3 and ADER-3 files retained for continued Report 2 development

## Current Layout

- `src/core`: configuration, state variables, fluxes, output, and solver dispatch
- `src/methods`: Godunov, MUSCL, WENO, and ADER method implementations
- `src/problems`: test problems
- `src/backends`: serial, OpenMP, and CUDA backend code
- `results`: selected Stage A--D outputs and figures
- `report_tex`: LaTeX source for the submitted Report 1 text

## Build In WSL

```bash
cd /mnt/d/code/Final_projects/report1_code_submission_2026-05-22
cmake -S . -B build
cmake --build build -j
./build/final_project
```

## Run Examples

Godunov:

```bash
./build/final_project --method=godunov --nx=200 --output=results/sod_1d_godunov.csv
```

MUSCL-Hancock:

```bash
./build/final_project --method=muscl --nx=200 --output=results/sod_1d_muscl.csv
```

WENO-2 + RK2:

```bash
./build/final_project --method=weno2 --nx=200 --output=results/sod_1d_weno2.csv
```

WENO-3 + RK3:

```bash
./build/final_project --method=weno3 --nx=200 --output=results/sod_1d_weno3.csv
```

Smooth advection + error norms:

```bash
./build/final_project --problem=smooth_advection --method=godunov --nx=200
./build/final_project --problem=smooth_advection --method=muscl --nx=200
```

Note:

- `smooth_advection` uses `final_time=0.3` by default
- This avoids error cancellation after a full-period return and is better suited for convergence-order checks.
- The default advection speed is `0.73`, and the smooth initial profile is phase shifted to avoid special alignment with the grid or period.

Gaussian advection:

```bash
./build/final_project --problem=gaussian_advection --method=godunov --nx=200
./build/final_project --problem=gaussian_advection --method=muscl --nx=200
./build/final_project --problem=gaussian_advection --method=weno2 --nx=200
```

Convergence study:

```bash
./build/final_project --study=convergence --problem=smooth_advection --method=godunov
./build/final_project --study=convergence --problem=smooth_advection --method=muscl
./build/final_project --study=convergence --problem=smooth_advection --method=weno2
./build/final_project --study=convergence --problem=smooth_advection --method=weno3
./build/final_project --study=convergence --problem=gaussian_advection --method=muscl
```

Read the summary csv:

```bash
python3 scripts/plot_convergence.py results/smooth_advection_1d_muscl_convergence.csv
```

Plot 2D fields:

```bash
python3 scripts/plot_field_2d.py results/quadrant_riemann_2d_godunov.csv
python3 scripts/plot_field_2d.py results/shock_bubble_2d_muscl.csv --field=rho
```

Render the main 6 Stage A 2D figures with unified style:

```bash
python3 scripts/render_stageA_2d_plots.py --dpi=300 --figsize 9 7
```

Compare two 2D fields with the same color scale:

```bash
python3 scripts/compare_fields_2d.py \
  results/quadrant_riemann_2d_godunov.csv \
  results/quadrant_riemann_2d_muscl.csv \
  --field rho
```

Compare serial vs openmp CSV results:

```bash
python3 scripts/compare_csv_results.py results/serial_muscl.csv results/openmp_muscl.csv
python3 scripts/compare_csv_results.py results/serial_2d_muscl.csv results/openmp_2d_muscl.csv
```

Run the Stage B discontinuous consistency checks:

```bash
python3 scripts/run_stageB_discontinuous_checks.py
```

CUDA scaffold:

```bash
cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON
cmake --build build-cuda -j
./build-cuda/final_project --backend=cuda
```

## Notes

- The main Report 1 results use Godunov, MUSCL-Hancock, and WENO-2.
- Preliminary WENO-3 and ADER-3 source files are included as continuing high-order development, but they are not required to reproduce the main Report 1 figures.
- `smooth_advection`, `gaussian_advection`, and `gaussian_advection_2d` use periodic boundaries and print density error norms when an exact reference is available.
- `study=convergence` automatically writes a summary CSV when it finishes
- 2D CSV outputs can be plotted directly with `scripts/plot_field_2d.py`.

