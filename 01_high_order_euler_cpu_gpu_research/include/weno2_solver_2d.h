#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// WENO-2 2D solver used for Report 1 accuracy and backend comparison studies.
SimulationResult2D run_weno2_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state);

}  // namespace fp
