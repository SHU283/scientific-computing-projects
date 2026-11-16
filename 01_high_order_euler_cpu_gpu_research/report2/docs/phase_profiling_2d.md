# Report 2 HLLC phase profiling

This diagnostic separates the smooth 2D solver runtime into setup, CFL reduction,
spatial reconstruction/flux work, state updates, admissibility checks, and transfers.
It is intended to identify where OpenMP and CUDA accelerate the solver and which
serial or synchronization-heavy phases limit the total speedup.

The profiler is opt-in. Set `FP_PHASE_PROFILE=1` to emit one machine-readable
`PHASE_PROFILE` line. Ordinary performance runs do not add the profiling
synchronizations and remain unchanged.

The canonical comparison uses `100 x 100` as the low-resolution case and `600 x
600` as the high-resolution case. The smaller case is large enough to avoid being
only a launch-overhead test, while the larger case matches the main Report 2
performance endpoint.

CUDA WENO-2 and WENO-3 use fused spatial kernels, so their reconstruction, HLLC
fluxes, and flux-divergence assembly are reported as `spatial_fused_s`. The CPU
paths expose x flux, y flux, and RHS assembly separately. The plotting script sums
these into the common `spatial operator` category. CUDA profiling deliberately
synchronizes between measured kernel groups; therefore its absolute profiled time
is diagnostic and should not replace the canonical no-I/O benchmark time.

The measured interval starts inside the solver. Command-line parsing, construction
of the initial condition, final error-norm evaluation, and CSV output are excluded.
`setup/transfer` includes host/device allocation, the initial host-to-device copy,
the final device-to-host copy, workspace allocation, and cleanup. The per-step copy
of the two reduced CUDA signal speeds is included in `CFL reduction`.
For cross-backend figures, state update and admissibility are combined as
`update/validity`: the CPU checks admissibility inside its update loop, while the
CUDA ADER-3 update kernel marks invalid cells and then copies one validity scalar
to the host. Their raw components remain separate in the CSV.

Run and plot with:

```bash
python3 scripts/run_report2_phase_profile_2d.py \
  --build-dir build-report2-hllc-cuda \
  --resolutions 100 600 --repeats 1 --resume
python3 scripts/plot_report2_phase_profile_2d.py \
  results/hllc/performance/phase_profile/phase_profile_2d_hllc_lovelace_a30.csv
python3 scripts/summarize_report2_phase_profile_2d.py \
  results/hllc/performance/phase_profile/phase_profile_2d_hllc_lovelace_a30.csv
```
