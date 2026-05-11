#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Experimental ADER-3 scalar-advection prototype rebuilt into Euler variables for output.
SimulationResult1D run_ader3_1d(const RuntimeConfig& config,
                                const std::vector<Conservative1D>& initial_state);

}  // namespace fp
