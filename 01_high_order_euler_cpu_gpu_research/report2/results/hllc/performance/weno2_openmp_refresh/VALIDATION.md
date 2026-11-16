# WENO-2 OpenMP refresh validation

Date: 2026-07-18

## Scope

This refresh follows the WENO-2 OpenMP alignment correction: the CFL maximum
reduction and both SSPRK2 state-update stages now use the OpenMP path.  Serial,
CUDA, WENO-3, and ADER-3 rows were not rerun or replaced by this campaign.

## Completeness checks

- Smooth Gaussian: 5 rows at `n = 25, 50, 100, 200, 600`.
- Discontinuous timing: 8 rows covering both cases at `n = 50, 100, 200, 600`.
- Thread scaling: 5 rows at `1, 2, 4, 8, 16` threads.
- Phase profile: 3 repeats at `100 x 100` and one diagnostic run at `600 x 600`.
- All rows use `method=weno2`, `backend=openmp`, the intended CFL, and zero
  phase-profile retries.
- Smooth error norms and step counts agree with the established numerical
  results, so the parallel correction did not change the solution.

## Timing quality

- Thread-scaling relative standard deviations are between 0.7% and 3.4%.
- The largest repeated discontinuous relative standard deviation is 8.9% for
  shock-bubble at `n=200`; its median remains consistent with the scaling trend.
- The two `n=600` discontinuous rows are deliberately single-run long-grid
  endpoints and must not be presented as repeated uncertainty estimates.
- Absolute end-to-end times differ from the older campaign because the runs
  occurred under a different shared-machine load.  The corrected rows are used
  for code-version consistency; repeat statistics and medians remain available
  for sensitivity checks.

## Phase-profile effect

For WENO-2 OpenMP at `600 x 600`, the corrected profile changes from the old
`25.86 s` total and `2.99 s` CFL cost to `19.51 s` total and `0.22 s` CFL cost.
The CFL fraction falls from about 11.6% to 1.1%, while the spatial-operator share
rises from 48.5% to 59.5%.  This is the expected signature of an active OpenMP
maximum reduction rather than a serial CFL scan.

## Canonical integration

`scripts/merge_weno2_openmp_refresh.py` replaced only matching WENO-2 OpenMP
rows in the canonical smooth, discontinuous, thread-scaling, and phase-profile
tables.  The pre-merge files are preserved under:

`results/report2_hllc/performance/archive/pre_weno2_openmp_refresh_20260718/`

