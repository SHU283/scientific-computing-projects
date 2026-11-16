#pragma once

#include <string>

#include <vector>

#include "config.h"
#include "euler1d.h"
#include "euler2d.h"

namespace fp {

// CUDA availability and runtime checks.
bool cuda_backend_compiled();
bool cuda_backend_probe();
std::string cuda_backend_summary();

// 1D CUDA solvers used by the Report 1 methods.
SimulationResult1D cuda_copy_1d_state(const RuntimeConfig& config,
                                      const std::vector<Conservative1D>& initial_state);
SimulationResult1D cuda_run_godunov_1d(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state);
SimulationResult1D cuda_run_muscl_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);
SimulationResult1D cuda_run_weno2_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state);

// 2D CUDA solvers. WENO-2, WENO-3, and ADER-3 perform the CFL reduction
// on the device, copying only the two directional maxima to the host.
SimulationResult2D cuda_run_godunov_2d(const RuntimeConfig& config,
                                       const Grid2D& grid,
                                       const Field2D& initial_state);
SimulationResult2D cuda_run_muscl_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);
SimulationResult2D cuda_run_weno2_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);
SimulationResult2D cuda_run_weno3_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);
SimulationResult2D cuda_run_ader3_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state);

}  // namespace fp
