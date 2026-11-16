# Formal Report 2 HLLC audit

Date: 2026-07-20

## Scope and authority

The formal report in `report/standalone_report2/` uses `code/hllc_solver/` and
`results/hllc/` as its authoritative implementation and result route.
HLL/Rusanov is not used for the abstract,
validation, backend-consistency evidence, main HLLC performance tables, or
official field figures. Recommendations and conclusions are based on HLLC;
they mention HLL only when summarising the explicitly matched comparison.

HLL numerical data appear only as a labelled historical comparator in the matched
accuracy--efficiency analysis. The two comparison figures and their table are
generated from processed strict matches under
`analysis/hll_vs_hllc/processed_results/`; they do not replace the HLLC
main dataset.

## Abstract values checked against canonical HLLC data

The smooth performance source is:

`results/hllc/performance/runtime_accuracy_2d_report2_hllc_lovelace_a30.csv`

At 600 x 600, the canonical HLLC values are:

| Method | L1 density error | Serial (s) | OpenMP (s) | CUDA (s) | OpenMP speedup | CUDA speedup |
|---|---:|---:|---:|---:|---:|---:|
| WENO-2 | 4.70263e-7 | 79.85748 | 19.56291 | 2.07785 | 4.0821 | 38.4327 |
| WENO-3 | 1.72624e-8 | 227.30518 | 37.13473 | 2.71128 | 6.1211 | 83.8369 |
| ADER-3 | 1.68102e-8 | 4118.48147 | 277.20461 | 10.50597 | 14.8572 | 392.0136 |

Therefore the abstract's rounded CUDA runtimes `2.08`, `2.71`, and `10.51 s`,
and speedups `83.84` and `392.01`, are HLLC values. The density-error reduction
factors relative to WENO-2 are `27.2420` for WENO-3 and `27.9749` for ADER-3,
reported as `27.2` and `28.0`.

The two-dimensional convergence source is:

`results/hllc/validation_2d_lovelace/report2_high_order_convergence_2d.csv`

The observed L1 orders over the finest 200-to-400 interval are `2.007007` for
WENO-2, `2.998085` for WENO-3, and `2.998018` for ADER-3. These are the values
rounded in the abstract and validation chapter.

The backend-consistency source is:

`results/hllc/backend_consistency_2d_lovelace/backend_consistency_2d_report2_lovelace_a30.csv`

It confirms zero stored OpenMP differences from serial, near-roundoff CUDA
differences for WENO, maxima below `4e-12` for smooth and quadrant ADER-3, and
the `2.0207295e-7` shock-bubble ADER-3 energy difference quoted in the report.

## Formal-report integration check

- Chapters on convergence, backend consistency, physical fields, discontinuous
  behaviour, phase profiling, OpenMP scaling, and GPU scaling use HLLC-labelled
  figures and HLLC canonical tables.
- The main 600 x 600 table in the accuracy--performance chapter matches the
  post-refresh HLLC CSV above.
- The only figures containing HLL data are
  `performance_hll_hllc_n600_tradeoff.png` and
  `performance_hll_hllc_cuda_tradeoff.png`, both inside the explicitly labelled
  matched HLL--HLLC comparison section.
- Discussion and conclusions treat HLLC as authoritative and refer to HLL only
  when summarising that matched comparison.

## Result

No HLL headline value was found in the abstract or the HLLC main tables. The
abstract now states the provenance explicitly so that the role of HLL cannot be
misread.

## CFL-aligned refresh (2026-07-21)

The formal validation and consistency sources were refreshed from the
consolidated CFL-aligned campaign. The official policy is:

| Problem class | WENO-2 | WENO-3 | ADER-3 |
|---|---:|---:|---:|
| Smooth 1D/2D | 0.40 | 0.40 | 0.25 |
| Discontinuous 1D/2D | 0.20 | 0.20 | 0.14 |

The canonical validation, backend-consistency, and high-resolution 1D paths now
contain the aligned results. Formal convergence, consistency, difference-map,
Sod, and 2D field figures were regenerated or refreshed from those paths. The
performance datasets were audited separately and already followed this policy.

Superseded pre-policy outputs were retained, not deleted, under
`_archive_completed_work_2026-07-19/superseded_cfl_results_2026-07-21/`.
They are excluded from the formal report and future submission material.
