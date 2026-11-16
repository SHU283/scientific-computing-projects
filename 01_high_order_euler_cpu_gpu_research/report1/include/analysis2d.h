#pragma once

#include <string>
#include <vector>

#include "euler2d.h"

namespace fp {

// Container for the three 1D error norms.
struct ErrorNorms2D {
    double l1_rho = 0.0;
    double l2_rho = 0.0;
    double linf_rho = 0.0;
};

// Compute 2D density error norms from a numerical field and a reference field.
ErrorNorms2D compute_density_error_2d(const Field2D& numerical, const Field2D& reference);
// Print 2D error norms to the terminal.
void print_error_norms_2d(const ErrorNorms2D& norms);
// Estimate the 2D convergence order.
double estimate_order_2d(double coarse_error, double fine_error);
// Print a 2D convergence table from grid resolutions and matching error norms.
void print_convergence_table_2d(const std::vector<int>& resolutions,
                                const std::vector<ErrorNorms2D>& norms);
// Write 2D convergence-study results to a CSV file.
void write_convergence_csv_2d(const std::string& path,
                              const std::vector<int>& resolutions,
                              const std::vector<ErrorNorms2D>& norms);

}  // namespace fp

