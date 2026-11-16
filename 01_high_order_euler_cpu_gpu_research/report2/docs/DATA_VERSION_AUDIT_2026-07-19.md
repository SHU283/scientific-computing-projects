# Report 2 HLLC data-version audit

Date: 2026-07-19

## Authoritative route

The official Report 2 implementation is under `code/hllc_solver/`, with its
curated numerical results under `results/hllc/`. HLL/Rusanov data remain a
comparison baseline only and are not
the source for Stage A validation, Stage B backend consistency, or the final
HLLC conclusions.

## July WENO-2 OpenMP refresh

The 18 July refresh corrected the WENO-2 OpenMP CFL reduction and both SSPRK2
state-update stages.  The following canonical performance tables contain the
replacement WENO-2 OpenMP rows:

- `results/hllc/performance/runtime_accuracy_2d_report2_hllc_lovelace_a30.csv`
- `results/hllc/performance/discontinuous/runtime_discontinuous_2d_report2_hllc_lovelace_a30_50_100_200.csv`
- `results/hllc/performance/discontinuous/runtime_discontinuous_2d_report2_hllc_lovelace_a30_n600.csv`
- `results/hllc/performance/discontinuous/runtime_discontinuous_2d_report2_hllc_lovelace_a30_50_100_200_600.csv`
- `results/hllc/performance/openmp_thread_scaling/openmp_thread_scaling_report2_hllc_lovelace_a30.csv`
- `results/hllc/performance/phase_profile/phase_profile_2d_hllc_lovelace_a30.csv`

The original files are preserved under
`results/hllc/performance/archive/pre_weno2_openmp_refresh_20260718/`.
They are historical records and must not be used for current figures or tables.

## Numerical-field verification after the correction

The Stage B Lovelace table was generated before the performance-only OpenMP
refresh.  To verify that its zero OpenMP differences remain valid for the
current source, the latest local executable was used to replay WENO-2 serial
and OpenMP at 100 x 100 for:

- Gaussian advection, CFL 0.4, final time 0.1: 48 steps;
- quadrant Riemann, CFL 0.2, final time 0.2: 213 steps;
- shock-bubble, CFL 0.2, final time 0.2: 427 steps.

For every case, the serial and OpenMP exported CSV files were byte-identical.
The audit outputs are stored separately under
`results/hllc/audit_latest_weno2_openmp_20260719/`; they do not replace
the canonical Lovelace backend-consistency CSV.

## Submission integration checked

- `report/standalone_report2/` uses the post-refresh smooth, discontinuous,
  thread-scaling, and phase-profile figures and values.
- `analysis/hll_vs_hllc/` contains the regenerated strict-match comparison
  tables and figures, including the accepted discontinuous runtime evidence.
- The historical HLL/Rusanov implementation is retained under
  `code/reference_hll/` and is used only as a labelled comparison baseline.
- Current HLLC headline values and figure provenance were re-audited on
  20 July 2026; see `FORMAL_REPORT_HLLC_AUDIT_2026-07-20.md`. No archived
  pre-refresh WENO-2 performance row is used in the formal report.

## Current canonical headline values

At 600 x 600 for smooth Gaussian advection, WENO-2 uses 79.86 s serial,
19.56 s OpenMP, and 2.08 s CUDA, corresponding to speedups of 4.08 and 38.43.
For the discontinuous 600 x 600 endpoint, the corrected WENO-2 OpenMP times are
95.18 s for quadrant Riemann and 173.33 s for shock-bubble.

The current WENO-2 16-thread OpenMP scaling speedup is 4.204.  At 600 x 600,
the corrected WENO-2 phase profile assigns 59.5 percent of OpenMP solver time
to the spatial operator and gives a spatial-operator speedup of 12.20 relative
to serial.
