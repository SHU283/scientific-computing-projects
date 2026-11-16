#pragma once

#include <string>
#include <vector>

#include "euler2d.h"

namespace fp {

struct ErrorNorms2D {
    double l1_rho = 0.0;
    double l2_rho = 0.0;
    double linf_rho = 0.0;
};

ErrorNorms2D compute_density_error_2d(const Field2D& numerical, const Field2D& reference);
void print_error_norms_2d(const ErrorNorms2D& norms);
double estimate_order_2d(double coarse_error, double fine_error);
void print_convergence_table_2d(const std::vector<int>& resolutions,
                                const std::vector<ErrorNorms2D>& norms);
void write_convergence_csv_2d(const std::string& path,
                              const std::vector<int>& resolutions,
                              const std::vector<ErrorNorms2D>& norms);

}  // namespace fp
