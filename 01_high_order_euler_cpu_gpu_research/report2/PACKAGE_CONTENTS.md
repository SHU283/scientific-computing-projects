# Package contents

## Implementation

`code/hllc_solver/` contains the complete source tree for the official HLLC
route. It includes the serial, OpenMP, and CUDA backends; WENO-2 with SSPRK2,
WENO-3 with SSPRK3, and ADER-3; one- and two-dimensional problems; tests; and
the scripts used for validation and analysis.

`code/reference_hll/` provides the HLL/Rusanov baseline for the matched
HLL--HLLC comparisons.

## Canonical HLLC results

| Directory | Purpose |
| --- | --- |
| `results/hllc/validation_1d/` | One-dimensional convergence and field validation. |
| `results/hllc/validation_2d_lovelace/` | Two-dimensional smooth validation on Lovelace. |
| `results/hllc/cfl_aligned_validation/` | Consolidated final-CFL validation and backend checks. |
| `results/hllc/extended_convergence_1600/` | Fine-grid convergence extension through 1600. |
| `results/hllc/backend_consistency_2d_lovelace/` | Serial/OpenMP/CUDA complete-field comparisons. |
| `results/hllc/high_resolution_fields/` | Accepted high-resolution physical fields and slices. |
| `results/hllc/discontinuous_evolution_2d/` | Time evolution for the two discontinuous 2D problems. |
| `results/hllc/performance/` | Smooth runtime and error campaign. |
| `results/hllc/performance/discontinuous/` | Discontinuous runtime and speedup campaign. |
| `results/hllc/performance/openmp_thread_scaling/` | OpenMP thread-count scaling. |
| `results/hllc/performance/phase_profile/` | Solver-phase timing and derived figures. |
| `results/hllc/performance/gpu_limit/` | Accepted CUDA-only scaling through 5000 x 5000. |

## ADER-3 CFL experiments

`results/supporting_experiments/ader3_cfl_selection/` contains the selected
ADER-3 CFL probe summaries, logs, and fields used to select the stability
settings. The main discontinuous experiments use CFL 0.14 for ADER-3.

## Figures and processed analysis

`figures/formal_report/` is a direct copy of the figure set used by the formal
Report 2 LaTeX project. `analysis/hll_vs_hllc/` contains the separately
processed matched comparison table, candidate figures, summary, rebuilding
script, and the eight HLL/HLLC input CSVs used by that script.

`report/standalone_report2/` contains `main.tex`, all section files, the crest,
and the same formal figure assets, so the standalone Report 2 document can be
compiled without referring back to the working repository. The combined
dissertation source is kept separately because it also packages Report 1.

## Documentation

- [Results](docs/results.md): convergence, backend agreement, performance, and regression checks.
- [HLLC implementation](docs/REPORT2_HLLC_ROUTE.md): flux construction and CFL selection.
- [Phase profiling](docs/phase_profiling_2d.md): timing categories and plotting commands.
