#pragma once

#include "config.h"
#include "euler2d.h"

namespace fp {

// Smooth 2D periodic transport benchmark and its exact shifted solution.
Field2D make_gaussian_advection_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);
Field2D make_gaussian_advection_2d_exact_solution(const RuntimeConfig& config,
                                                  const Grid2D& grid,
                                                  double time);

// Discontinuous 2D benchmark problems for qualitative field comparison.
Field2D make_quadrant_riemann_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);
Field2D make_shock_bubble_2d_initial_condition(const RuntimeConfig& config, const Grid2D& grid);

}  // namespace fp
