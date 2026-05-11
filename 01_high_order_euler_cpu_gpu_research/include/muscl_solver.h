#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Second-order MUSCL-Hancock 1D solver with minmod limiting and HLL fluxes.
SimulationResult1D run_muscl_hancock_1d(const RuntimeConfig& config,
                                        const std::vector<Conservative1D>& initial_state);

}  // namespace fp
