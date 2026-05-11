#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// Second-order 2D MUSCL solver with directional limited reconstruction and SSP-RK2.
SimulationResult2D run_muscl_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state);

}  // namespace fp
