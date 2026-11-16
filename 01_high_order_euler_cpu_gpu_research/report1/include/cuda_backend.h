#pragma once

#include <string>

#include <vector>

#include "config.h"
#include "euler1d.h"
#include "euler2d.h"

namespace fp {

// CUDA environment checks.
// Return whether the CUDA backend was compiled into the program.
bool cuda_backend_compiled();
// Probe whether the current runtime environment can actually use the CUDA backend.
bool cuda_backend_probe();
// Return a short text description of the CUDA backend status.
std::string cuda_backend_summary();

// 1D CUDA solver interface.
// Copy a 1D initial state to the CUDA backend and return the result object.
SimulationResult1D cuda_copy_1d_state(const RuntimeConfig& config,
                                      const std::vector<Conservative1D>& initial_state);
// Run the 1D Godunov method on the CUDA backend.
SimulationResult1D cuda_run_godunov_1d(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state);
// Run the 1D MUSCL method on the CUDA backend.
SimulationResult1D cuda_run_muscl_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);
// Run the 1D WENO-2 method on the CUDA backend.
SimulationResult1D cuda_run_weno2_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);

// 2D CUDA solver interface.
// Run the 2D Godunov method on the CUDA backend.
SimulationResult2D cuda_run_godunov_2d(const RuntimeConfig& config,
                                       const Grid2D& grid,
                                       const Field2D& initial_state);
// Run the 2D MUSCL method on the CUDA backend.
SimulationResult2D cuda_run_muscl_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);
// Run the 2D WENO-2 method on the CUDA backend.
SimulationResult2D cuda_run_weno2_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);

}  // namespace fp

