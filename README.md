# Selected Scientific Computing Projects

## 中文快速说明

本仓库包含三个与科学计算、CFD 数值仿真、高性能计算和 GPU/并行计算相关的代表性项目，主要用于展示数值方法实现、并行计算实践和实验结果整理能力。

项目主线是：基于 C++ 实现可压缩 Euler 方程有限体积求解器，并围绕不同数值格式、CPU/OpenMP/MPI/CUDA 后端、运行时间、加速比、结果一致性和可视化后处理进行实验分析。

建议快速查看顺序：

1. `01_high_order_euler_cpu_gpu_research`
   - 当前毕业论文/研究项目，最能体现科学计算、数值方法、C++、OpenMP、CUDA 和性能分析能力。
   - 重点看 `README.md`、`src/methods/`、`src/backends/`、`scripts/` 和 `results/`。

2. `02_shock_bubble_cuda_cpu_gpu_solver`
   - 重点体现 CUDA 版本实现、CPU/GPU 结果对比、运行时间对比和二维激波-气泡算例。
   - 侧重展示 GPU/仿真方向的实现与结果对比。

3. `03_shock_bubble_openmp_mpi_solver`
   - 重点体现从串行程序扩展到 OpenMP 与 MPI 并行版本，包括域分解、halo exchange、强/弱扩展测试和并行效率分析。
   - 侧重展示 HPC/并行计算方向的实现与性能分析。

说明：完整中间数据和大规模输出文件较大，因此这里只保留代表性结果、图像、CSV 和脚本；源码、构建方式和复现实验入口均保留。

This repository contains three selected scientific-computing / CFD projects.

## Projects

1. `01_high_order_euler_cpu_gpu_research`
   - Ongoing MPhil Scientific Computing project.
   - Finite-volume methods for compressible Euler equations.
   - Serial CPU, OpenMP, and CUDA backend comparison.
   - Includes representative accuracy, runtime, backend-consistency, and field-visualization results.

2. `02_shock_bubble_cuda_cpu_gpu_solver`
   - CUDA/CPU comparison project for a 2D shock-bubble benchmark.
   - Godunov-HLL and MUSCL-Hancock finite-volume methods.
   - Includes representative CPU/GPU output fields, difference plots, runtime plots, and selected raw CSV samples.

3. `03_shock_bubble_openmp_mpi_solver`
   - OpenMP/MPI parallel finite-volume solver for a 2D shock-bubble benchmark.
   - Includes domain decomposition, halo exchange, strong/weak scaling scripts, figures, logs, and selected raw CSV samples.

## Notes

- Generated build folders, IDE caches, and Python bytecode caches are excluded.
- Large full time-series output folders are not fully included; representative raw data and plots are provided instead.
- Each project has its own `README.md`, `PROJECT_NOTE.md`, and result notes.
