# Finite-Volume Euler Solver — Report 1 and Report 2

Curated submission packages for Shuli Xu's two project reports on finite-volume
Euler solvers, numerical accuracy, and serial/OpenMP/CUDA performance.

## Repository guide

| Package | Contents | Start here |
| --- | --- | --- |
| `report1/` | Godunov, MUSCL-Hancock, and WENO-2 source; experiment scripts; selected Stage A–D data and figures | [Report 1 submission guide](report1/README_SUBMISSION.md) |
| `report2/` | Official HLLC solver; HLL/Rusanov comparison source; canonical data; figures; analysis; standalone Report 2 LaTeX source | [Report 2 guide](report2/README.md) and [package contents](report2/PACKAGE_CONTENTS.md) |

Report 2 extends the study to WENO-3 and ADER-3 and includes convergence,
backend consistency, HLL–HLLC comparisons, phase profiling, OpenMP scaling,
and the accepted CUDA scaling campaign through 5000 × 5000.

## Build

Use CMake 3.20 or later and a C++17 compiler from Linux or WSL.

Report 1:

```bash
cd report1
cmake -S . -B build
cmake --build build --parallel
./build/final_project --problem=gaussian_advection_2d --method=muscl --nx=100 --ny=100
```

Report 2:

```bash
cd report2/code/hllc_solver
cmake -S . -B build-report2 -DFP_ENABLE_OPENMP=ON -DFP_ENABLE_CUDA=OFF
cmake --build build-report2 --parallel
python3 scripts/run_report2_local_smoke.py --build-dir build-report2
```

CUDA builds require a compatible CUDA toolkit and GPU; configure with
`-DFP_ENABLE_CUDA=ON`. Plotting dependencies are listed in
[`report2/requirements.txt`](report2/requirements.txt).

## Reports and remote experiments

Use [Report 1's build and experiment guide](report1/README_SUBMISSION.md) to
run its solver. The standalone Report 2 LaTeX source is under
`report2/report/standalone_report2/`.

The Report 2 remote workflow scripts use Lovelace-specific workspace and
account paths. Set those paths for your own machine before running a remote
campaign. The analysis and plotting scripts run against the packaged data.

## 中文说明

项目包含 Report 1 和 Report 2 的数值方法实现、并行计算代码、实验数据和图表。
Report 1 从 `report1/README_SUBMISSION.md` 开始阅读；Report 2 从
`report2/README.md` 开始阅读。Report 2 采用 HLLC 通量，HLL/Rusanov 用于对照实验。
