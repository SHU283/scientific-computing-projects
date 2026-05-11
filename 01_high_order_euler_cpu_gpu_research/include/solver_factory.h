#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Select the requested 1D numerical method from RuntimeConfig::method.
SimulationResult1D run_selected_solver(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state);

}  // namespace fp
