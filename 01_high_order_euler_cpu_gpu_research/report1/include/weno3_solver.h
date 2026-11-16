#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Run the 1D WENO-3 solver with RK3 time integration.
SimulationResult1D run_weno3_rk3_1d(const RuntimeConfig& config,
                                    const std::vector<Conservative1D>& initial_state);

}  // namespace fp

