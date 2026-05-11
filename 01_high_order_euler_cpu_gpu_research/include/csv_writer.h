#pragma once

#include <string>

#include "euler1d.h"

namespace fp {

// Write a 1D cell-centered solution as x,rho,u,p,energy for post-processing.
void write_result_csv(const std::string& path, const SimulationResult1D& result, double gamma);

}  // namespace fp
