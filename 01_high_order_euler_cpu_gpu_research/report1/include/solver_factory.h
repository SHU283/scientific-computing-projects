#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Parse, validate, or apply runtime configuration settings.
// One-dimensional problem setup or update step.
SimulationResult1D run_selected_solver(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state);

}  // namespace fp

