#pragma once

#include <string>

#include "euler2d.h"

namespace fp {

void write_result_csv_2d(const std::string& path, const SimulationResult2D& result, double gamma);

}  // namespace fp
