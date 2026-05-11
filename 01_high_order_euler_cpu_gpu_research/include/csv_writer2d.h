#pragma once

#include <string>

#include "euler2d.h"

namespace fp {

// Write a 2D cell-centered field as i,j,rho,u,v,p,energy for plotting scripts.
void write_result_csv_2d(const std::string& path, const SimulationResult2D& result, double gamma);

}  // namespace fp
