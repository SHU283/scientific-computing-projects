# Selected Scientific Computing Projects

## 中文快速说明

本仓库包含三个与科学计算、CFD 数值仿真、高性能计算和 GPU/并行计算相关的代表性项目，主要用于展示数值方法实现、并行计算实践和实验结果整理能力。

项目主线是：基于 C++ 实现可压缩 Euler 方程有限体积求解器，并围绕不同数值格式、CPU/OpenMP/MPI/CUDA 后端、运行时间、加速比、结果一致性和可视化后处理进行实验分析。

建议快速查看顺序：

1. [01_high_order_euler_cpu_gpu_research](01_high_order_euler_cpu_gpu_research/)
   - 毕业论文/研究项目，包含 Report 1 和 Report 2 的完整整理版提交包，展示数值方法、C++、OpenMP、CUDA 和性能分析。
   - [Report 1](01_high_order_euler_cpu_gpu_research/report1/)：Godunov、MUSCL-Hancock、WENO-2，及 Stage A–D 验证、后端一致性和性能结果。
   - [Report 2](01_high_order_euler_cpu_gpu_research/report2/)：HLLC 主实现、HLL/Rusanov 对照、WENO-2/WENO-3/ADER-3、正式实验数据、图表和独立报告 LaTeX 源码。
   - 从[构建与复现说明](01_high_order_euler_cpu_gpu_research/README.md)开始阅读；包内附逐文件 SHA-256 校验清单。

2. `02_shock_bubble_cuda_cpu_gpu_solver`
   - 重点体现 CUDA 版本实现、CPU/GPU 结果对比、运行时间对比和二维激波-气泡算例。
   - 侧重展示 GPU/仿真方向的实现与结果对比。

3. `03_shock_bubble_openmp_mpi_solver`
   - 重点体现从串行程序扩展到 OpenMP 与 MPI 并行版本，包括域分解、halo exchange、强/弱扩展测试和并行效率分析。
   - 侧重展示 HPC/并行计算方向的实现与性能分析。

说明：项目 01 保留 Report 1 和 Report 2 已整理提交包的全部文件，包含正式采用的数据、图表、代码和复现说明。项目 02、03 保留代表性结果、图像、CSV 和脚本。构建产物、缓存及未纳入提交包的中间运行结果不在展示范围内。

This repository contains three selected scientific-computing / CFD projects.

## Projects

1. [01_high_order_euler_cpu_gpu_research](01_high_order_euler_cpu_gpu_research/)
   - MPhil Scientific Computing research project: complete curated Report 1 and Report 2 submission packages.
   - Finite-volume methods for compressible Euler equations.
   - Serial CPU, OpenMP, and CUDA backend comparison.
   - [Report 1](01_high_order_euler_cpu_gpu_research/report1/): Godunov, MUSCL-Hancock, WENO-2, and Stage A–D validation and performance evidence.
   - [Report 2](01_high_order_euler_cpu_gpu_research/report2/): the official HLLC route, the HLL/Rusanov comparison, WENO-2/WENO-3/ADER-3, canonical data, figures, and standalone report source.
   - Includes the package build guides and a file-by-file SHA-256 manifest.

2. `02_shock_bubble_cuda_cpu_gpu_solver`
   - CUDA/CPU comparison project for a 2D shock-bubble benchmark.
   - Godunov-HLL and MUSCL-Hancock finite-volume methods.
   - Includes representative CPU/GPU output fields, difference plots, runtime plots, and selected raw CSV samples.

3. `03_shock_bubble_openmp_mpi_solver`
   - OpenMP/MPI parallel finite-volume solver for a 2D shock-bubble benchmark.
   - Includes domain decomposition, halo exchange, strong/weak scaling scripts, figures, logs, and selected raw CSV samples.

## Notes

- Generated build folders, IDE caches, and Python bytecode caches are excluded.
- Project 01 preserves all 1,123 original files from the two curated submission packages. Additional development archives and superseded campaigns remain outside their scope.
- Projects 02 and 03 contain representative raw data and plots rather than every large time-series output.
- Each project has its own `README.md`; project-specific submission guides or project/result notes provide further details.
