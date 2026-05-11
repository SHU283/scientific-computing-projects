#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// First-order 2D finite-volume solver using Rusanov fluxes on Cartesian faces.
SimulationResult2D run_godunov_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state);

}  // namespace fp
