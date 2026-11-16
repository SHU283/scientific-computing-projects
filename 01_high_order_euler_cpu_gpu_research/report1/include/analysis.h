#pragma once

#include <string>
#include <vector>

#include "euler1d.h"

namespace fp {

// Container for the three 1D error norms.
struct ErrorNorms1D {
    double l1_rho = 0.0;
    double l2_rho = 0.0;
    double linf_rho = 0.0;
};

// Compute 1D density error norms from a numerical solution and a reference solution.
ErrorNorms1D compute_density_error(const SimulationResult1D& numerical,
                                   const SimulationResult1D& reference);
// Print 1D error norms to the terminal.
void print_error_norms(const ErrorNorms1D& norms);
// Estimate the 1D convergence order.
double estimate_order(double coarse_error, double fine_error);
// Print a 1D convergence table from grid resolutions and matching error norms.
void print_convergence_table(const std::vector<int>& resolutions,
                             const std::vector<ErrorNorms1D>& norms);
// Write 1D convergence-study results to a CSV file.
void write_convergence_csv(const std::string& path,
                           const std::vector<int>& resolutions,
                           const std::vector<ErrorNorms1D>& norms);

}  // namespace fp

