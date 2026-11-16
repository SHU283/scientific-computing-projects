#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// MUSCL reconstruction step using limited slopes.
SimulationResult2D run_muscl_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state);

}  // namespace fp

