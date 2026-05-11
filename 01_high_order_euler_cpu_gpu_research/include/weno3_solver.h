#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Experimental WENO-3 + SSP-RK3 extension for 1D smooth-advection studies.
SimulationResult1D run_weno3_rk3_1d(const RuntimeConfig& config,
                                    const std::vector<Conservative1D>& initial_state);

}  // namespace fp
