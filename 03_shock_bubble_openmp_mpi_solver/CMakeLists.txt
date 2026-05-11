cmake_minimum_required(VERSION 3.16)
project(hrsc_euler LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(serial
    src/main_serial.cpp
    src/field.cpp
    src/init_shockbubble.cpp
    src/bc.cpp
    src/dt.cpp
    src/hrsc_muscl_hancock.cpp
    src/output.cpp
    src/checks.cpp
)

target_include_directories(serial PRIVATE include)
target_compile_options(serial PRIVATE -O3 -march=native -Wall -Wextra)
target_link_libraries(serial PRIVATE MPI::MPI_CXX)

add_executable(omp
    src/main_omp.cpp
    src/field.cpp
    src/init_shockbubble.cpp
    src/bc.cpp
    src/dt.cpp
    src/hrsc_muscl_hancock.cpp
    src/output.cpp
    src/checks.cpp
)

target_include_directories(omp PRIVATE include)
target_compile_options(omp PRIVATE -O3 -march=native -Wall -Wextra -fopenmp)
target_link_options(omp PRIVATE -fopenmp)
target_link_libraries(omp PRIVATE OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(mpi
  src/main_mpi.cpp
  src/mpi_decomp.cpp
  src/mpi_output.cpp
  src/hrsc_muscl_hancock_mpi.cpp
  src/hrsc_muscl_hancock.cpp
  src/field.cpp
  src/init_shockbubble.cpp
  src/bc.cpp
  src/dt.cpp
  src/output.cpp
  src/checks.cpp
)

target_include_directories(mpi PRIVATE include)
target_link_libraries(mpi PRIVATE MPI::MPI_CXX)
