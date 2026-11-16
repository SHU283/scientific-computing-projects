# Report 2 submission evidence audit

Audit date: 2026-08-04; final technical verification: 2026-08-05

## Verdict

The curated package contains the active standalone Report 2 source, every
figure it references, the canonical CSV or raw-field evidence behind its
quantitative tables, and the scripts/data needed for the principal processed
comparisons. It also retains useful supporting experiment records without
mixing them into the authoritative HLLC result tree.

The formal route is HLLC. HLL/Rusanov is retained only as the matched baseline
for the accuracy-efficiency comparison. The official GPU-limit sequence ends
at 5000 x 5000.

## Report-to-evidence mapping

| Report content | Packaged evidence |
| --- | --- |
| 1D smooth/Gaussian convergence and orders | `results/hllc/extended_convergence_1600/validation_1d_openmp/` |
| 1D smooth, Gaussian, and Sod profiles | `results/hllc/cfl_aligned_validation/high_resolution_fields_1d/` and `results/hllc/validation_1d/` |
| 2D Gaussian convergence through 1600 | `results/hllc/extended_convergence_1600/validation_2d_cuda/` |
| 600 x 600 Gaussian, quadrant-Riemann, and shock-bubble fields/slices | `results/hllc/high_resolution_fields/` |
| Quadrant-Riemann and shock-bubble time evolution | `results/hllc/discontinuous_evolution_2d/` |
| Serial/OpenMP/CUDA field consistency and difference maps | `results/hllc/backend_consistency_2d_lovelace/` and `results/hllc/cfl_aligned_validation/` |
| Positivity/admissibility statements for resolved fields | Canonical exported fields in the high-resolution and consistency directories |
| ADER-3 discontinuous CFL-selection table/discussion | `results/supporting_experiments/ader3_cfl_selection/` |
| Smooth runtime, error, speedup, and accuracy-runtime plots | `results/hllc/performance/runtime_accuracy_2d_report2_hllc_lovelace_a30.csv` |
| Discontinuous runtime and speedup | `results/hllc/performance/discontinuous/` |
| HLL-HLLC matched comparison | `analysis/hll_vs_hllc/source_data/`, `processed_results/`, and `scripts/build_hll_vs_hllc_analysis.py` |
| Solver-phase breakdown and phase speedup | `results/hllc/performance/phase_profile/` |
| OpenMP thread-count scaling | `results/hllc/performance/openmp_thread_scaling/` |
| CUDA-only scaling and error refinement | `results/hllc/performance/gpu_limit/` |

The evidence-matrix and final recommendation tables are synthesis tables in
the LaTeX source. They combine the mapped experiments above rather than being
independent raw measurements.

## Figure integrity

- All 30 graphics referenced by the active standalone Report 2 source resolve
  inside `report/standalone_report2/`.
- The formal figure directory contains the 29 report PNGs; the thirtieth
  graphic is the University crest PDF.
- All formal PNG files were decoded successfully during this audit.
- `figures/formal_report/` mirrors the report figure assets for direct review.

## Supporting and reusable research material

Beyond the minimum report figures, the package includes complete-field CSVs,
intermediate discontinuous states, convergence extensions, repeated timing
rows, phase-level profiles, OpenMP thread sweeps, accepted large-grid CUDA
measurements, and the selected CFL probes. These records can support further
questions about spatial localisation of backend differences, contact
resolution, scaling exponents, phase bottlenecks, and sensitivity to CFL.

The HLL-HLLC comparison is self-contained: its local `source_data/` directory
contains the eight immutable input CSVs used by the packaged rebuilding script.
Processed tables are saved separately and do not overwrite those inputs.

## Canonical boundaries and exclusions

The following are deliberately not treated as submission evidence:

- exploratory 6200/7000-grid timings with non-monotone shared-node behaviour;
- failed or interrupted reruns;
- compiler/build output and caches;
- transfer archives and duplicate packaging files;
- broad superseded-CFL campaigns, except the small labelled selection subset.

The standalone Report 2 LaTeX project is included. The combined Report 1--2
dissertation project remains outside this technical package because it also
contains the previously submitted Report 1 material.

## Reproduction notes

The C++17/OpenMP/CUDA source is in `code/hllc_solver/`. The accepted raw and
processed results are immutable package inputs. Plotting requirements are
listed in the root `requirements.txt`. A local LaTeX PDF was not rebuilt during
this audit; the included source and assets are arranged as an Overleaf-ready
standalone project.

## Final technical verification

The following checks were repeated on 2026-08-05 after the submission package
was curated and its comments were edited:

- Both solver trees built successfully in serial-only, OpenMP, and CUDA
  configurations, giving six successful CMake builds in total.
- A 68-case executable smoke matrix completed successfully. It covered HLLC
  and the retained HLL/Rusanov reference, serial/OpenMP/CUDA backends,
  WENO-2/WENO-3/ADER-3, the smooth two-dimensional Gaussian problem, both
  discontinuous two-dimensional problems, and representative one-dimensional
  paths at the formal CFL values.
- The packaged backend-alignment audit reran all nine HLLC
  problem--method combinations at 100 x 100. Serial and 16-thread OpenMP used
  identical step counts, and every exported field agreed exactly in the stored
  decimal output (`max_delta = 0`).
- All 72 Python files passed an AST syntax check, all 16 shell scripts passed
  `bash -n`, and all 35 PowerShell scripts passed the PowerShell parser.
- All 254 CSV/TSV files were structurally readable with consistent row widths.
  All 123 PNG files and both PDF files had valid file signatures; no zero-byte
  files, build directories, bytecode caches, or editor caches remained.
- The HLL--HLLC analysis was regenerated in an isolated temporary directory.
  Its nine processed CSV files were byte-identical to the packaged versions,
  and all 15 candidate comparison figures were reproduced.
- The phase-profile audit, formal CFL-policy audit, and formal report-asset
  audit all passed. The active LaTeX source resolves 30 graphics and has closed
  cross-reference and citation sets.

No TeX engine is installed in the local verification environment, so a fresh
PDF could not be compiled there. This is the only unexecuted build step. The
active source graph, bibliography, labels, references, and graphic paths were
instead checked statically, and the project remains ready for its normal
two-pass `pdflatex` build on Overleaf or another TeX installation.
