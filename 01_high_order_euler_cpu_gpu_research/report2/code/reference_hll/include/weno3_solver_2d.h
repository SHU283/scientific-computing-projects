#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

SimulationResult2D run_weno3_rk3_2d(const RuntimeConfig& config,
                                    const Grid2D& grid,
                                    const Field2D& initial_state);

}  // namespace fp
