# HRSC Euler Solver with OpenMP and MPI

## 中文快速说明

这是一个二维可压缩 Euler 方程激波-气泡相互作用项目，重点是将串行有限体积求解器扩展为 OpenMP 共享内存并行版本和 MPI 分布式内存并行版本。

项目主要能力点：

- 基于 C++ 实现二维 Euler 方程有限体积求解器。
- 使用 MUSCL-Hancock、minmod limiter、HLL Riemann solver、Strang splitting 等数值方法。
- OpenMP 版本用于多线程循环并行。
- MPI 版本包含域分解、halo exchange、rank 间边界通信、结果汇总输出。
- 使用强扩展、弱扩展、运行时间、加速比和并行效率图评估并行效果。

推荐查看顺序：

1. `src/`：串行、OpenMP、MPI 版本的核心实现。
2. MPI decomposition / halo-exchange 相关代码：重点体现分布式并行能力。
3. `figs/`：代表性场图与强/弱扩展性能图。
4. `scripts/`：画图和 scaling 结果整理脚本。

说明：本项目重点不在复杂物理模型本身，而在从串行数值程序扩展到 OpenMP/MPI 并行求解器，并完成基础性能分析。

## Overview

This project implements a two-dimensional finite-volume solver for the compressible Euler equations. The main benchmark is the classical shock-bubble interaction problem, where a Mach 1.22 shock wave in air interacts with a helium bubble.

The code extends a serial solver to OpenMP shared-memory parallelism and MPI distributed-memory parallelism. It is mainly used to study solver correctness, parallelization structure, and scaling behavior.

## Numerical Methods

- MUSCL-Hancock second-order reconstruction.
- Minmod slope limiter.
- HLL approximate Riemann solver.
- Strang dimensional splitting.
- Transmissive boundary conditions.
- CFL-based time-step control.

## Repository Structure

- `include/`: headers for field storage, parameters, boundary conditions, time-step calculation, solver routines, MPI decomposition, output, and checks.
- `src/`: serial, OpenMP, and MPI implementations.
- `scripts/`: Python scripts for plotting solution fields and scaling results.
- `figs/`: selected generated figures.
- `output/`: generated CSV output directory, excluded from the repository.
- `CMakeLists.txt`: build configuration for serial, OpenMP, and MPI targets.

## Build

```bash
mkdir -p build
cd build
cmake ..
make -j
```

## Notes for Review

- Large CSV outputs and build artifacts are intentionally excluded.
- The recommended review entry points are the solver routines in `src/` and the MPI decomposition / halo-exchange code.
- The project emphasizes the transition from serial implementation to OpenMP and MPI versions, including runtime, speedup, and parallel-efficiency analysis.
