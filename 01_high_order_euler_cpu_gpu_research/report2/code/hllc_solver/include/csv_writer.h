#pragma once

#include <string>

#include "euler1d.h"

namespace fp {

void write_result_csv(const std::string& path, const SimulationResult1D& result, double gamma);

}  // namespace fp
