#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// Run the 2D WENO-2 solver with SSP-RK2 time integration.
SimulationResult2D run_weno2_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state);

}  // namespace fp

