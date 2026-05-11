cmake_minimum_required(VERSION 3.20)

project(FinalProjectEulerSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FP_ENABLE_OPENMP "Enable OpenMP backend" ON)
option(FP_ENABLE_CUDA "Enable CUDA backend" OFF)

if(FP_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
    else()
        message(WARNING "FP_ENABLE_CUDA=ON but no CUDA compiler was found. Building without CUDA.")
        set(FP_ENABLE_CUDA OFF)
    endif()
endif()

add_executable(final_project
    src/main.cpp
    src/core/analysis.cpp
    src/core/analysis2d.cpp
    src/core/config.cpp
    src/core/euler1d.cpp
    src/core/euler2d.cpp
    src/core/csv_writer.cpp
    src/core/csv_writer2d.cpp
    src/core/solver_factory.cpp
    src/methods/godunov_solver.cpp
    src/methods/godunov_solver_2d.cpp
    src/methods/muscl_solver.cpp
    src/methods/muscl_solver_2d.cpp
    src/methods/ader3_solver.cpp
    src/methods/weno2_solver.cpp
    src/methods/weno2_solver_2d.cpp
    src/methods/weno3_solver.cpp
    src/problems/gaussian_advection_2d.cpp
    src/problems/sod_shock_tube.cpp
)

target_include_directories(final_project
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

if(FP_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(final_project PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(final_project PRIVATE FP_HAS_OPENMP=1)
    endif()
endif()

if(FP_ENABLE_CUDA)
    target_sources(final_project PRIVATE src/backends/cuda_backend.cu)
    set_target_properties(final_project PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
    target_compile_definitions(final_project PRIVATE FP_HAS_CUDA=1)
else()
    target_sources(final_project PRIVATE src/backends/cuda_backend_stub.cpp)
endif()

message(STATUS "OpenMP enabled: ${FP_ENABLE_OPENMP}")
message(STATUS "CUDA enabled: ${FP_ENABLE_CUDA}")
