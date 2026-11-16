# Official Report 2 HLLC solver

This is the authoritative Report 2 implementation. It supports serial,
OpenMP, and CUDA execution for WENO-2/SSPRK2, WENO-3/SSPRK3, and ADER-3, with
the full HLLC Euler interface flux used by the Report 2 methods.

Formal CFL values are 0.40 for smooth WENO-2/WENO-3, 0.25 for smooth ADER-3,
0.20 for discontinuous WENO-2/WENO-3, and 0.14 for discontinuous ADER-3.

## Build

```bash
cmake -S . -B build-report2 -DFP_ENABLE_OPENMP=ON -DFP_ENABLE_CUDA=OFF
cmake --build build-report2 --parallel
```

Enable CUDA with `-DFP_ENABLE_CUDA=ON` on a machine with a compatible CUDA
toolkit. `scripts/run_report2_local_smoke.py` provides a compact regression
suite after building.

Generated build directories and result files are intentionally absent from
this source directory. Canonical outputs are stored at
`../../results/hllc/` in the submission package.

