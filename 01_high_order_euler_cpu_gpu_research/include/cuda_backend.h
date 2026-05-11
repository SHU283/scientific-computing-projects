#pragma once

#include <string>

#include <vector>

#include "config.h"
#include "euler1d.h"
#include "euler2d.h"

namespace fp {

// Build/runtime status helpers used by the command-line driver before CUDA runs.
bool cuda_backend_compiled();
bool cuda_backend_probe();
std::string cuda_backend_summary();

// Diagnostic copy-only path for validating basic GPU memory transfer.
SimulationResult1D cuda_copy_1d_state(const RuntimeConfig& config,
                                      const std::vector<Conservative1D>& initial_state);

// CUDA 1D solver variants used for CPU/GPU consistency and runtime comparisons.
SimulationResult1D cuda_run_godunov_1d(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state);
SimulationResult1D cuda_run_muscl_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);
SimulationResult1D cuda_run_weno2_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);

// CUDA 2D solver variants for Report 1 smooth and discontinuous benchmark cases.
SimulationResult2D cuda_run_godunov_2d(const RuntimeConfig& config,
                                       const Grid2D& grid,
                                       const Field2D& initial_state);
SimulationResult2D cuda_run_muscl_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);
SimulationResult2D cuda_run_weno2_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);

}  // namespace fp
