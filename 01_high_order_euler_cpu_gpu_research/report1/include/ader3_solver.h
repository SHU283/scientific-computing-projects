#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Run the 1D ADER-3 solver.
SimulationResult1D run_ader3_1d(const RuntimeConfig& config,
                                const std::vector<Conservative1D>& initial_state);

}  // namespace fp

