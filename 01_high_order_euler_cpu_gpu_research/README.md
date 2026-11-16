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

## Snapshot integrity and scope

The 450 files under `report1/` and 673 files under `report2/` preserve the
submission package contents byte for byte. They were checked against the
original Report 1 `.tar.gz` and Report 2 `.zip` archives before import.
[`SUBMISSION_MANIFEST.csv`](SUBMISSION_MANIFEST.csv) records each packaged
file's size and SHA-256 checksum. Git line-ending conversion is disabled to
preserve these archived bytes across platforms.

The original package documentation is retained. In particular, the historical
Report 1 `README.md` mentions `report_tex/` and an old absolute checkout path;
the submitted archive does not contain `report_tex/`. Use the build commands
above and `report1/README_SUBMISSION.md` for this repository. The standalone
Report 2 LaTeX source is available under `report2/report/standalone_report2/`.

The Report 2 remote workflow scripts retain their original local and Lovelace
paths as provenance. Adjust those paths for another machine or account before
rerunning remote campaigns. The package-local analysis remains self-contained.

The original compressed archives, separate signed declarations, presentation
materials, and wider development workspace are outside this repository's
submission-package scope. Existing technical validation is documented in
[`report2/docs/SUBMISSION_EVIDENCE_AUDIT_2026-08-04.md`](report2/docs/SUBMISSION_EVIDENCE_AUDIT_2026-08-04.md).

## 中文说明

本仓库分别保存 Report 1 和 Report 2 已整理好的提交包。代码、数据、图表和
包内文档保持原始内容，并附文件校验清单。Report 1 从
`report1/README_SUBMISSION.md` 开始阅读；Report 2 从 `report2/README.md`
开始阅读。Report 2 的 HLLC 路线是正式实现，HLL/Rusanov 仅作为对照。
