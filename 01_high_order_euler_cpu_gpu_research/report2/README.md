# Report 2 submission package

This directory is the curated reproducibility package for Report 2. It contains
the source code, canonical data, report figures, analysis scripts, and concise
provenance documentation needed to inspect or reproduce the submitted results.

## Authoritative route

HLLC is the authoritative interface flux for Report 2. The main implementation
is under `code/hllc_solver/`, and the main datasets are under `results/hllc/`.
The older HLL/Rusanov source under `code/reference_hll/` is retained only for
the matched HLL--HLLC accuracy-efficiency comparison.

The formal CFL policy is:

| Problem class | WENO-2 | WENO-3 | ADER-3 |
| --- | ---: | ---: | ---: |
| Smooth | 0.40 | 0.40 | 0.25 |
| Discontinuous | 0.20 | 0.20 | 0.14 |

For a given experiment row, serial, OpenMP, and CUDA use the same method,
problem, grid, final time, boundary condition, and CFL value.

## Directory guide

- `code/hllc_solver/`: official Report 2 C++17/OpenMP/CUDA implementation.
- `code/reference_hll/`: older HLL/Rusanov implementation used as a comparator.
- `results/hllc/`: canonical validation, field, consistency, performance,
  profiling, OpenMP scaling, and accepted GPU-limit data.
- `results/supporting_experiments/`: narrowly selected non-canonical evidence
  retained to document experimental decisions, currently the ADER-3 CFL study.
- `figures/formal_report/`: the figures used by the formal Report 2 source.
- `report/standalone_report2/`: self-contained standalone Report 2 LaTeX source
  retained for review and word-count checking.
- `analysis/hll_vs_hllc/`: processed matched tables, comparison figures, and
  a self-contained script plus its immutable HLL/HLLC input CSVs.
- `docs/`: CFL, provenance, figure-inventory, and result-version records.
- `workflows/`: retained remote workflow scripts for the official campaigns.
See `PACKAGE_CONTENTS.md` for the detailed scope.

## Build and smoke test

From `code/hllc_solver/`, a CPU/OpenMP build can be produced with:

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

The plotting and comparison scripts use the Python packages listed in
`requirements.txt`.

The scripts under `workflows/` preserve the original Lovelace collection
environment, including the Windows workspace root and remote account paths.
They are provenance records for the submitted campaigns. To rerun them from a
different checkout or account, update their explicit local and remote roots;
the package-local analysis and audit scripts do not depend on those paths.

Implementation comments in both solver trees and the retained workflow scripts
are written in English. Comments focus on numerical choices, backend behaviour,
stability safeguards, and reproducibility rather than restating individual
statements.

## Curation policy

Original raw records have not been overwritten. This package contains the
accepted results through the official 5000 x 5000 CUDA endpoint. Superseded
CFL campaigns are excluded from the canonical result tree. A small, clearly
labelled subset of ADER-3 CFL probes is retained under
`results/supporting_experiments/ader3_cfl_selection/` solely to support the
reported CFL-selection discussion. Exploratory 6200/7000-grid runs, failed
reruns, compiler output, build directories, caches, and compressed transfer
archives are intentionally excluded. The formal source-to-result mapping is
recorded in `docs/FORMAL_RESULT_PROVENANCE_2026-07-21.md`; the submission-wide
coverage check is in `docs/SUBMISSION_EVIDENCE_AUDIT_2026-08-04.md`.

The combined Report 1--Report 2 dissertation project is intentionally kept
outside this technical package because it also contains the already submitted
Report 1 material. The standalone Report 2 source is included here.
