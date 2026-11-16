#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

Field2D make_gaussian_advection_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);
Field2D make_gaussian_advection_2d_exact_solution(const RuntimeConfig& config,
                                                  const Grid2D& grid,
                                                  double time);
Field2D make_quadrant_riemann_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);
Field2D make_shock_bubble_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);

}  // namespace fp

