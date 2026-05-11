# 2D Compressible Euler Solver for Shock-Bubble Interaction

## 中文快速说明

这是一个二维可压缩 Euler 方程激波-气泡相互作用算例项目，重点是比较 CPU 与 CUDA GPU 版本在结果一致性和运行时间上的差异。

项目主要能力点：

- 使用 C++/CUDA 实现二维有限体积求解器。
- 围绕 shock-bubble benchmark 组织数值实验。
- 实现 Godunov-HLL 与 MUSCL-Hancock 等格式，并进行 CFL 时间步控制和 ghost-cell 边界处理。
- 使用 Python 脚本生成密度/压力/速度场图像、CPU/GPU 差异图和性能图。

推荐查看顺序：

1. `src/main/main.cpp`：CPU 主流程与算例入口。
2. `src/gpu/gpu_main.cu`：CUDA 版本主流程与 GPU kernel 调用。
3. `project_figs/`：代表性场图、差异图和性能图。
4. `tools/`：结果后处理与画图脚本。

说明：本项目更偏 CUDA/GPU 对比和二维 CFD benchmark 实现，用于展示 GPU 仿真开发、结果一致性检查和基础性能分析流程。

## Overview

This project implements a two-dimensional finite-volume solver for the compressible Euler equations. The main benchmark is the shock-bubble interaction problem, where a Mach 1.22 shock wave propagates through air and interacts with a helium bubble.

The project contains both CPU and CUDA implementations and is used to compare numerical accuracy, runtime behavior, and CPU/GPU result consistency.

## Numerical Methods

- Godunov finite-volume method with HLL flux.
- MUSCL-Hancock second-order reconstruction.
- Strang splitting for the second-order update.
- Ghost-cell boundary treatment.
- CFL-based time-step calculation.

## Repository Structure

- `include/`: header files.
- `src/`: CPU, GPU, initialization, output, and driver code.
- `tools/`: Python scripts for plotting and post-processing.
- `project_figs/`: selected generated figures.
- `data/`: generated CSV output directory, excluded from the repository.

## Build

```bash
mkdir -p build
cd build
cmake ..
make -j
```

Generated fields are written to a local `data/` directory when the solver is run.

## Notes for Review

- Large CSV outputs and build artifacts are intentionally excluded.
- The recommended review entry points are `src/main/main.cpp`, `src/gpu/gpu_main.cu`, and the plotting scripts in `tools/`.
- The project focuses on finite-volume implementation, CPU/GPU comparison, and basic performance analysis rather than production CFD software packaging.
