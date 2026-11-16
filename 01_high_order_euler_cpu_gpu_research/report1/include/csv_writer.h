#pragma once

#include <string>

#include "euler1d.h"

namespace fp {

// Write a 1D result to a CSV file.
void write_result_csv(const std::string& path, const SimulationResult1D& result, double gamma);

}  // namespace fp

