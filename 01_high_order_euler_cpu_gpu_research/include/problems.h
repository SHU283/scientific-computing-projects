#pragma once

#include <vector>

#include "config.h"
#include "euler1d.h"

namespace fp {

// Dispatch to the configured 1D benchmark initial condition.
std::vector<Conservative1D> make_initial_condition(const RuntimeConfig& config);

// 1D benchmark families used in convergence and shock-capturing tests.
std::vector<Conservative1D> make_sod_initial_condition(const RuntimeConfig& config);
std::vector<Conservative1D> make_smooth_advection_initial_condition(const RuntimeConfig& config);
std::vector<Conservative1D> make_gaussian_advection_initial_condition(const RuntimeConfig& config);

// Exact solution for smooth periodic advection benchmarks.
SimulationResult1D make_exact_solution(const RuntimeConfig& config, double time);

}  // namespace fp
