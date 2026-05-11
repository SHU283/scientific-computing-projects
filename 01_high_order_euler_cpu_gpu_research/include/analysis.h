#pragma once

#include <string>
#include <vector>

#include "euler1d.h"

namespace fp {

// Error norms for density-based 1D convergence studies.
struct ErrorNorms1D {
    double l1_rho = 0.0;
    double l2_rho = 0.0;
    double linf_rho = 0.0;
};

// Compare a numerical solution against a reference solution on the same grid.
ErrorNorms1D compute_density_error(const SimulationResult1D& numerical,
                                   const SimulationResult1D& reference);
void print_error_norms(const ErrorNorms1D& norms);

// Utility functions for convergence-rate reporting.
double estimate_order(double coarse_error, double fine_error);
void print_convergence_table(const std::vector<int>& resolutions,
                             const std::vector<ErrorNorms1D>& norms);
void write_convergence_csv(const std::string& path,
                           const std::vector<int>& resolutions,
                           const std::vector<ErrorNorms1D>& norms);

}  // namespace fp
