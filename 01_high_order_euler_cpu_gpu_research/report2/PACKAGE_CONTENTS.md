# Package contents

## Submitted implementation

`code/hllc_solver/` contains the complete source tree for the official HLLC
route. It includes the serial, OpenMP, and CUDA backends; WENO-2 with SSPRK2,
WENO-3 with SSPRK3, and ADER-3; one- and two-dimensional problems; tests; and
the scripts used for validation and analysis.

`code/reference_hll/` is a preserved comparison baseline. It is not the
authoritative Report 2 solver and should only be used to reproduce the matched
HLL--HLLC comparisons in Section 3 of the analysis.

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

## Supporting decision evidence

`results/supporting_experiments/ader3_cfl_selection/` contains the selected
ADER-3 CFL probe summaries, logs, and fields needed to audit the stability
choice described in the report. These files are supporting experimental
history, not an alternative canonical dataset. The formal discontinuous
results use CFL 0.14 for ADER-3.

## Figures and processed analysis

`figures/formal_report/` is a direct copy of the figure set used by the formal
Report 2 LaTeX project. `analysis/hll_vs_hllc/` contains the separately
processed matched comparison table, candidate figures, summary, rebuilding
script, and the eight immutable HLL/HLLC source CSVs used by that script. The
analysis is therefore reproducible inside the package and does not silently
replace either route's raw records.

`report/standalone_report2/` contains `main.tex`, all section files, the crest,
and the same formal figure assets, so the standalone Report 2 document can be
compiled without referring back to the working repository. The combined
dissertation source is kept separately because it also packages Report 1.

## Documentation

- `docs/FORMAL_RESULT_PROVENANCE_2026-07-21.md`: authoritative provenance map.
- `docs/CFL_AUDIT_2026-07-20.md`: formal CFL audit.
- `docs/DATA_VERSION_AUDIT_2026-07-19.md`: data-version review.
- `docs/FORMAL_REPORT_HLLC_AUDIT_2026-07-20.md`: formal-report route audit.
- `docs/REPORT2_FIGURE_INVENTORY.md`: figure and evidence inventory.
- `docs/phase_profiling_2d.md`: interpretation of phase-profile categories.
- `docs/SUBMISSION_EVIDENCE_AUDIT_2026-08-04.md`: package-wide traceability and
  completeness audit.

## Deliberate exclusions

The package does not include build products, `.git` metadata, Python caches,
transfer archives, failed diagnostics, or the later 6200/7000-grid exploratory
runs. Superseded CFL data are also excluded except for the explicitly labelled
selection evidence described above. Those other files remain in the working
repository but are not part of the formal submission evidence.
