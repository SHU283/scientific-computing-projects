#include "solver_factory.h"

#include <stdexcept>

#include "godunov_solver.h"
#include "muscl_solver.h"
#include "ader3_solver.h"
#include "weno2_solver.h"
#include "weno3_solver.h"

namespace fp {

SimulationResult1D run_selected_solver(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state) {
    if (config.method == "godunov") {
        return run_godunov_1d(config, initial_state);
    }

    if (config.method == "muscl") {
        return run_muscl_hancock_1d(config, initial_state);
    }

    if (config.method == "weno2") {
        return run_weno2_rk2_1d(config, initial_state);
    }

    if (config.method == "weno3") {
        return run_weno3_rk3_1d(config, initial_state);
    }

    if (config.method == "ader3") {
        return run_ader3_1d(config, initial_state);
    }

    throw std::runtime_error("Unknown method: " + config.method);
}

}  // namespace fp
