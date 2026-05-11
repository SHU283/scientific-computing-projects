cmake_minimum_required(VERSION 3.18)
project(cfd_cuda_new LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")

set(CMAKE_CUDA_ARCHITECTURES 86)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(cpu
  src/main/main.cpp
  src/cpu/bc.cpp
  src/cpu/dt.cpp
  src/cpu/muscl.cpp
  src/cpu/riemann_hll.cpp
  src/cpu/step.cpp
  src/physics/init.cpp
  src/io/output.cpp
)

add_executable(gpu
  src/gpu/gpu_main.cu
  src/gpu/gpu_kernels.cu
  src/gpu/gpu_muscl.cu
  src/gpu/gpu_mem.cu
  src/gpu/gpu_bc.cu
  src/gpu/gpu_dt.cu
  src/physics/init.cpp
  src/io/output.cpp
)
