#include "analysis.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fp {

ErrorNorms1D compute_density_error(const SimulationResult1D& numerical,
                                   const SimulationResult1D& reference) {
    if (numerical.cells.size() != reference.cells.size()) {
        throw std::runtime_error("Numerical and reference solutions must have the same size.");
    }

    ErrorNorms1D norms;
    const double inv_n = 1.0 / static_cast<double>(numerical.cells.size());

    // Accumulate mean absolute, RMS, and maximum density errors.
    for (std::size_t i = 0; i < numerical.cells.size(); ++i) {
        const double err = std::abs(numerical.cells[i].rho - reference.cells[i].rho);
        norms.l1_rho += err;
        norms.l2_rho += err * err;
        norms.linf_rho = std::max(norms.linf_rho, err);
    }

    norms.l1_rho *= inv_n;
    norms.l2_rho = std::sqrt(norms.l2_rho * inv_n);
    return norms;
}

void print_error_norms(const ErrorNorms1D& norms) {
    std::cout << "density error L1   : " << norms.l1_rho << "\n";
    std::cout << "density error L2   : " << norms.l2_rho << "\n";
    std::cout << "density error Linf : " << norms.linf_rho << "\n";
}

double estimate_order(double coarse_error, double fine_error) {
    if (coarse_error <= 0.0 || fine_error <= 0.0) {
        return 0.0;
    }
    return std::log2(coarse_error / fine_error);
}

// Convergence-study bookkeeping and reporting step.
void print_convergence_table(const std::vector<int>& resolutions,
                             const std::vector<ErrorNorms1D>& norms) {
    if (resolutions.size() != norms.size()) {
        throw std::runtime_error("Convergence table input sizes do not match.");
    }

    // Print run information or diagnostics for the user.
    std::cout << "\nConvergence table (density)\n";
    // Print run information or diagnostics for the user.
    std::cout << "nx\tL1\t\torder(L1)\tL2\t\torder(L2)\tLinf\t\torder(Linf)\n";

    // First row has no observed order because there is no coarser reference.
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        std::cout << resolutions[i] << "\t"
                  << std::scientific << std::setprecision(6) << norms[i].l1_rho << "\t";

        if (i == 0) {
            std::cout << "-\t\t";
        } else {
            std::cout << estimate_order(norms[i - 1].l1_rho, norms[i].l1_rho) << "\t\t";
        }

        std::cout << norms[i].l2_rho << "\t";
        if (i == 0) {
            std::cout << "-\t\t";
        } else {
            std::cout << estimate_order(norms[i - 1].l2_rho, norms[i].l2_rho) << "\t\t";
        }

        std::cout << norms[i].linf_rho << "\t";
        if (i == 0) {
            std::cout << "-";
        } else {
            std::cout << estimate_order(norms[i - 1].linf_rho, norms[i].linf_rho);
        }
        std::cout << "\n";
    }

    std::cout << std::defaultfloat;
}

// Convergence-study bookkeeping and reporting step.
void write_convergence_csv(const std::string& path,
                           const std::vector<int>& resolutions,
                           const std::vector<ErrorNorms1D>& norms) {
    if (resolutions.size() != norms.size()) {
        throw std::runtime_error("Convergence csv input sizes do not match.");
    }

    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open convergence csv file: " + path);
    }

    out << "nx,l1_rho,order_l1,l2_rho,order_l2,linf_rho,order_linf\n";

    // Write one line per resolution, including observed orders relative to the
    // preceding coarser grid.
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        const double order_l1 = (i == 0) ? 0.0 : estimate_order(norms[i - 1].l1_rho, norms[i].l1_rho);
        const double order_l2 = (i == 0) ? 0.0 : estimate_order(norms[i - 1].l2_rho, norms[i].l2_rho);
        const double order_linf = (i == 0) ? 0.0 : estimate_order(norms[i - 1].linf_rho, norms[i].linf_rho);

        out << resolutions[i] << ","
            << norms[i].l1_rho << ","
            << order_l1 << ","
            << norms[i].l2_rho << ","
            << order_l2 << ","
            << norms[i].linf_rho << ","
            << order_linf << "\n";
    }
}

}  // namespace fp

