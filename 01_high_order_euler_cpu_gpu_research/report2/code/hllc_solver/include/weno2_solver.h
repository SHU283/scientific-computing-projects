#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

SimulationResult1D run_weno2_rk2_1d(const RuntimeConfig& config,
                                    const std::vector<Conservative1D>& initial_state);

}  // namespace fp
