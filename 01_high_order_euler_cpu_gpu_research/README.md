# High-Order Finite-Volume Euler Solver

## 中文快速说明

这是当前 MPhil Scientific Computing 毕业论文/研究项目的提交版，主题是面向二维可压缩 Euler 方程的高阶有限体积方法与多后端性能对比。

项目核心目标是：在统一 C++ 框架下实现 Godunov-HLL、MUSCL-Hancock、WENO-2 等数值格式，并比较 serial CPU、OpenMP 和 CUDA 后端在精度、运行时间、加速比和 CPU/GPU 结果一致性上的表现。

项目技术要点：

- C++ 数值程序开发：守恒量/原始量转换、Euler flux、HLL/Rusanov 数值通量、CFL 时间步控制。
- CFD/数值方法：有限体积法、Godunov、MUSCL、WENO、Sod shock tube、quadrant Riemann、shock-bubble 等算例。
- 并行与 GPU：OpenMP 后端、CUDA kernel、CPU/GPU 一致性检查、运行时间与加速比分析。
- 工程化后处理：CSV 输出、Python 画图、误差表、runtime 表、代表性结果整理。

推荐查看顺序：

1. `src/core/`：Euler 方程基础数据结构、通量和配置。
2. `src/methods/`：Godunov、MUSCL、WENO-2 等核心数值格式实现。
3. `src/backends/`：CUDA/OpenMP 相关后端实现。
4. `scripts/`：代表性实验、后端一致性检查和性能图表生成脚本。
5. `results/`：已筛选的代表性结果、图像和 CSV 表格。

说明：WENO-3 和 ADER-3 属于更长期论文计划中的实验性扩展，当前主线重点是 Godunov-HLL、MUSCL-Hancock 和 WENO-2。

## Overview

This repository contains an ongoing MPhil Scientific Computing project on finite-volume methods for the compressible Euler equations. The project compares different numerical schemes and execution backends on one-dimensional and two-dimensional test problems.

The current focus is on Godunov-HLL, MUSCL-Hancock, and WENO-2 methods, with serial CPU, OpenMP, and CUDA backends used for accuracy, runtime, and CPU/GPU consistency studies.

## Implemented Scope

- One-dimensional and two-dimensional compressible Euler test framework.
- Godunov-HLL, MUSCL-Hancock, and WENO-2 methods.
- Serial CPU, OpenMP, and CUDA backend structure.
- Smooth advection, Sod shock tube, quadrant Riemann, and shock-bubble style test cases.
- CSV output, Python post-processing, convergence summaries, and runtime comparison scripts.

Some higher-order extensions such as WENO-3 and ADER-3 are part of the broader thesis plan and may be incomplete or experimental.

## Repository Structure

- `include/`: public headers.
- `src/core/`: configuration, conserved variables, fluxes, output, and solver dispatch.
- `src/methods/`: finite-volume numerical methods.
- `src/problems/`: initial conditions and benchmark problems.
- `src/backends/`: serial, OpenMP, and CUDA backend code.
- `scripts/`: compact representative experiment drivers, backend checks, plotting, and runtime analysis.
- `tests/`: lightweight tests and validation utilities.
- `examples/`: selected figures and compact CSV summaries for quick review.
- `results/`: curated representative outputs and figures.
- `data/`, `build/`, and `build-cuda/`: local/generated artifacts excluded from the repository.

## Build

```bash
cmake -S . -B build
cmake --build build -j
./build/final_project
```

CUDA build, if CUDA is available:

```bash
cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON
cmake --build build-cuda -j
./build-cuda/final_project --backend=cuda
```

## Example Runs

```bash
./build/final_project --method=godunov --nx=200 --output=results/sod_1d_godunov.csv
./build/final_project --method=muscl --nx=200 --output=results/sod_1d_muscl.csv
./build/final_project --method=weno2 --nx=200 --output=results/sod_1d_weno2.csv
```

Representative script examples:

```bash
python3 scripts/run_report1_stageA_representatives.py
python3 scripts/run_report1_stageB_checks.py
python3 scripts/render_stageA_2d_plots.py
```

## Example Results

A small curated result set is included under `examples/`:

- method-comparison figures for smooth and discontinuous 2D cases,
- runtime and speedup summaries for serial, OpenMP, and CUDA backends,
- OpenMP thread-scaling and backend-consistency summaries,
- compact CSV tables used to generate the example plots.

Full field-output CSV files are intentionally excluded because they are large and can be regenerated from the source code.

## Notes for Review

- This is a research codebase rather than a production library.
- Only representative generated results are included; full intermediate output folders are intentionally excluded from this curated repository version.
- Recommended review entry points are `src/core/`, `src/methods/`, `src/backends/`, and the analysis scripts in `scripts/`.
