#include "analysis2d.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fp {

ErrorNorms2D compute_density_error_2d(const Field2D& numerical, const Field2D& reference) {
    if (numerical.nx != reference.nx || numerical.ny != reference.ny) {
        throw std::runtime_error("2D numerical and reference fields must have the same size.");
    }

    ErrorNorms2D norms;
    const double inv_n = 1.0 / static_cast<double>(numerical.cells.size());
    for (std::size_t k = 0; k < numerical.cells.size(); ++k) {
        const double err = std::abs(numerical.cells[k].rho - reference.cells[k].rho);
        norms.l1_rho += err;
        norms.l2_rho += err * err;
        norms.linf_rho = std::max(norms.linf_rho, err);
    }

    norms.l1_rho *= inv_n;
    norms.l2_rho = std::sqrt(norms.l2_rho * inv_n);
    return norms;
}

void print_error_norms_2d(const ErrorNorms2D& norms) {
    std::cout << "2D density error L1   : " << norms.l1_rho << "\n";
    std::cout << "2D density error L2   : " << norms.l2_rho << "\n";
    std::cout << "2D density error Linf : " << norms.linf_rho << "\n";
}

double estimate_order_2d(double coarse_error, double fine_error) {
    if (coarse_error <= 0.0 || fine_error <= 0.0) {
        return 0.0;
    }
    return std::log2(coarse_error / fine_error);
}

void print_convergence_table_2d(const std::vector<int>& resolutions,
                                const std::vector<ErrorNorms2D>& norms) {
    if (resolutions.size() != norms.size()) {
        throw std::runtime_error("2D convergence table input sizes do not match.");
    }

    std::cout << "\nConvergence table 2D (density)\n";
    std::cout << "n\tL1\t\torder(L1)\tL2\t\torder(L2)\tLinf\t\torder(Linf)\n";

    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        std::cout << resolutions[i] << "\t"
                  << std::scientific << std::setprecision(6) << norms[i].l1_rho << "\t";

        if (i == 0) {
            std::cout << "-\t\t";
        } else {
            std::cout << estimate_order_2d(norms[i - 1].l1_rho, norms[i].l1_rho) << "\t\t";
        }

        std::cout << norms[i].l2_rho << "\t";
        if (i == 0) {
            std::cout << "-\t\t";
        } else {
            std::cout << estimate_order_2d(norms[i - 1].l2_rho, norms[i].l2_rho) << "\t\t";
        }

        std::cout << norms[i].linf_rho << "\t";
        if (i == 0) {
            std::cout << "-";
        } else {
            std::cout << estimate_order_2d(norms[i - 1].linf_rho, norms[i].linf_rho);
        }
        std::cout << "\n";
    }

    std::cout << std::defaultfloat;
}

void write_convergence_csv_2d(const std::string& path,
                              const std::vector<int>& resolutions,
                              const std::vector<ErrorNorms2D>& norms) {
    if (resolutions.size() != norms.size()) {
        throw std::runtime_error("2D convergence csv input sizes do not match.");
    }

    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open 2D convergence csv file: " + path);
    }

    out << "n,l1_rho,order_l1,l2_rho,order_l2,linf_rho,order_linf\n";
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        const double order_l1 = (i == 0) ? 0.0 : estimate_order_2d(norms[i - 1].l1_rho, norms[i].l1_rho);
        const double order_l2 = (i == 0) ? 0.0 : estimate_order_2d(norms[i - 1].l2_rho, norms[i].l2_rho);
        const double order_linf = (i == 0) ? 0.0 : estimate_order_2d(norms[i - 1].linf_rho, norms[i].linf_rho);

        out << resolutions[i] << ','
            << norms[i].l1_rho << ','
            << order_l1 << ','
            << norms[i].l2_rho << ','
            << order_l2 << ','
            << norms[i].linf_rho << ','
            << order_linf << '\n';
    }
}

}  // namespace fp
