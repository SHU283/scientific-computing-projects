# Formal Report Result Provenance (2026-07-21)

## Official policy

| Problem class | WENO-2 | WENO-3 | ADER-3 |
|---|---:|---:|---:|
| Smooth 1D/2D | 0.40 | 0.40 | 0.25 |
| Discontinuous 1D/2D | 0.20 | 0.20 | 0.14 |

Serial, OpenMP, and CUDA use the same CFL for a matched method/problem row.

## Canonical data

- Convergence: `results/hllc/validation_1d/` and
  `results/hllc/validation_2d_lovelace/`.
- Backend consistency: `results/hllc/backend_consistency_2d_lovelace/`.
- High-resolution fields: `results/hllc/high_resolution_fields/`.
- Runtime, speedup, thread scaling, phase profiling, and accepted GPU-limit
  evidence: `results/hllc/performance/`.
- Strict HLL--HLLC comparison: `analysis/hll_vs_hllc/`.
- Immutable consolidated campaign:
  `results/hllc/cfl_aligned_validation/`.

## Formal figure provenance

The formal report copies are under `report/standalone_report2/figures/`. Convergence,
backend summary, backend difference maps, high-resolution fields, and all
performance figures are byte-identical to their canonical sources. The two
midline slices are byte-identical to the HLLC weekly-record figures generated
from the accepted 600 x 600 field exports.

Run both audits after any result or report refresh:

```text
python code/hllc_solver/scripts/audit_report2_cfl_policy.py
python code/hllc_solver/scripts/audit_formal_report_assets.py
```

The first checks CFL constants in active scripts and all self-describing
canonical CSV rows. The second checks that every formal figure exists, verifies
the official CFL statement in the methodology chapter, and confirms that the
formal and curated figure files have identical contents.

## Historical material

The ADER-3 CFL probe folders are retained as selection evidence only. They are
not canonical result sources. Superseded pre-policy validation and field outputs
were moved without deletion to:

```text
_archive_completed_work_2026-07-19/superseded_cfl_results_2026-07-21/
```

Nothing in that archive should enter the formal report or submission package.
