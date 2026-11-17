# Report 2: High-order Euler solvers on CPU and GPU

C++17 implementations of WENO-2, WENO-3, and ADER-3 with serial, OpenMP,
and CUDA backends. The project studies convergence, backend agreement,
flow fields, runtime scaling, and phase-level bottlenecks.

## Methods and experiment settings

The main solver uses HLLC interface fluxes and lives in `code/hllc_solver/`.
`code/reference_hll/` provides the HLL/Rusanov baseline for matched
accuracy-efficiency comparisons.

| Problem class | WENO-2 | WENO-3 | ADER-3 |
| --- | ---: | ---: | ---: |
| Smooth | 0.40 | 0.40 | 0.25 |
| Discontinuous | 0.20 | 0.20 | 0.14 |

For each backend comparison, serial, OpenMP, and CUDA use the same method,
problem, grid, final time, boundary condition, and CFL.

## Directory guide

- `code/hllc_solver/`: C++17/OpenMP/CUDA implementation and experiment scripts.
- `code/reference_hll/`: HLL/Rusanov comparison solver.
- `results/hllc/`: convergence, fields, backend comparisons, runtime,
  phase profiling, OpenMP scaling, and CUDA scaling through 5000 × 5000.
- `results/supporting_experiments/ader3_cfl_selection/`: probes used to select
  the discontinuous ADER-3 CFL.
- `figures/formal_report/`: report figures.
- `report/standalone_report2/`: standalone Report 2 LaTeX source.
- `analysis/hll_vs_hllc/`: matched comparison tables, figures, input CSVs,
  and the rebuilding script.
- `docs/`: [results](docs/results.md), [HLLC implementation](docs/REPORT2_HLLC_ROUTE.md),
  and [phase profiling](docs/phase_profiling_2d.md).
- `workflows/`: scripts for remote experiment campaigns.

See [package contents](PACKAGE_CONTENTS.md) for the result directories.

## Build and run

From `code/hllc_solver/`:

```bash
cmake -S . -B build-report2 \
  -DFP_ENABLE_OPENMP=ON \
  -DFP_ENABLE_CUDA=OFF
cmake --build build-report2 --parallel
python3 scripts/run_report2_local_smoke.py --build-dir build-report2
```

On a CUDA-capable machine with a compatible toolkit:

```bash
cmake -S . -B build-report2-cuda \
  -DFP_ENABLE_OPENMP=ON \
  -DFP_ENABLE_CUDA=ON
cmake --build build-report2-cuda --parallel
```

Install plotting dependencies from `requirements.txt`.
Remote scripts under `workflows/` use the original Lovelace workspace and
account paths; update their local and remote roots for your environment.
