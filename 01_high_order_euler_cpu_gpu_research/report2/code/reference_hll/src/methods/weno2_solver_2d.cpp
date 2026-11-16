#include "weno2_solver_2d.h"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef FP_HAS_OPENMP
#include <omp.h>
#endif

#ifndef FP_HAS_OPENMP
#define FP_HAS_OPENMP 0
#endif

namespace fp {

namespace {

constexpr double kWenoEps = 1.0e-6;

int wrap_index(int idx, int size) {
    while (idx < 0) {
        idx += size;
    }
    while (idx >= size) {
        idx -= size;
    }
    return idx;
}

const Conservative2D& sample_cell(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    if (config.boundary == "periodic") {
        return field(wrap_index(i, field.nx), wrap_index(j, field.ny));
    }

    i = std::max(0, std::min(i, field.nx - 1));
    j = std::max(0, std::min(j, field.ny - 1));
    return field(i, j);
}

double nonlinear_weight_z(double beta, double tau, double linear_weight) {
    const double ratio = tau / (kWenoEps + beta);
    return linear_weight * (1.0 + ratio * ratio);
}

double mapped_weight(double omega, double linear_weight) {
    const double d = linear_weight;
    const double numerator = omega * (d + d * d - 3.0 * d * omega + omega * omega);
    const double denominator = d * d + omega * (1.0 - 2.0 * d);
    if (std::abs(denominator) < 1.0e-14) {
        return omega;
    }
    return numerator / denominator;
}

double weno2_left_component(double um1, double u0, double up1) {
    const double beta0 = (u0 - um1) * (u0 - um1);
    const double beta1 = (up1 - u0) * (up1 - u0);
    const double tau = std::abs(beta0 - beta1);

    double alpha0 = nonlinear_weight_z(beta0, tau, 0.5);
    double alpha1 = nonlinear_weight_z(beta1, tau, 0.5);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = mapped_weight(omega0, 0.5);
    omega1 = mapped_weight(omega1, 0.5);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    const double q0 = -0.5 * um1 + 1.5 * u0;
    const double q1 = 0.5 * u0 + 0.5 * up1;
    return omega0 * q0 + omega1 * q1;
}

double weno2_right_component(double u0, double up1, double up2) {
    const double beta0 = (up1 - u0) * (up1 - u0);
    const double beta1 = (up2 - up1) * (up2 - up1);
    const double tau = std::abs(beta0 - beta1);

    double alpha0 = nonlinear_weight_z(beta0, tau, 0.5);
    double alpha1 = nonlinear_weight_z(beta1, tau, 0.5);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = mapped_weight(omega0, 0.5);
    omega1 = mapped_weight(omega1, 0.5);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    const double q0 = 0.5 * u0 + 0.5 * up1;
    const double q1 = 1.5 * up1 - 0.5 * up2;
    return omega0 * q0 + omega1 * q1;
}

Conservative2D reconstruct_left_state_x(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    const Conservative2D& um1 = sample_cell(field, i - 1, j, config);
    const Conservative2D& u0 = sample_cell(field, i, j, config);
    const Conservative2D& up1 = sample_cell(field, i + 1, j, config);
    return {
        weno2_left_component(um1.rho, u0.rho, up1.rho),
        weno2_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
        weno2_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
        weno2_left_component(um1.energy, u0.energy, up1.energy),
    };
}

Conservative2D reconstruct_right_state_x(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    const Conservative2D& u0 = sample_cell(field, i, j, config);
    const Conservative2D& up1 = sample_cell(field, i + 1, j, config);
    const Conservative2D& up2 = sample_cell(field, i + 2, j, config);
    return {
        weno2_right_component(u0.rho, up1.rho, up2.rho),
        weno2_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
        weno2_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
        weno2_right_component(u0.energy, up1.energy, up2.energy),
    };
}

Conservative2D reconstruct_left_state_y(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    const Conservative2D& um1 = sample_cell(field, i, j - 1, config);
    const Conservative2D& u0 = sample_cell(field, i, j, config);
    const Conservative2D& up1 = sample_cell(field, i, j + 1, config);
    return {
        weno2_left_component(um1.rho, u0.rho, up1.rho),
        weno2_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
        weno2_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
        weno2_left_component(um1.energy, u0.energy, up1.energy),
    };
}

Conservative2D reconstruct_right_state_y(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    const Conservative2D& u0 = sample_cell(field, i, j, config);
    const Conservative2D& up1 = sample_cell(field, i, j + 1, config);
    const Conservative2D& up2 = sample_cell(field, i, j + 2, config);
    return {
        weno2_right_component(u0.rho, up1.rho, up2.rho),
        weno2_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
        weno2_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
        weno2_right_component(u0.energy, up1.energy, up2.energy),
    };
}

double compute_time_step(const Field2D& field, const Grid2D& grid, const RuntimeConfig& config) {
    double max_speed = 1.0e-12;
    for (int j = 0; j < field.ny; ++j) {
        for (int i = 0; i < field.nx; ++i) {
            max_speed = std::max(max_speed, max_signal_speed_x(field(i, j), config.gamma));
            max_speed = std::max(max_speed, max_signal_speed_y(field(i, j), config.gamma));
        }
    }
    return config.cfl * std::min(dx(grid), dy(grid)) / max_speed;
}

Field2D add_scaled(const Field2D& base, const Field2D& delta, double scale) {
    Field2D result = base;
    for (std::size_t k = 0; k < base.cells.size(); ++k) {
        result.cells[k] = base.cells[k] + scale * delta.cells[k];
    }
    return result;
}

// Assemble x- and y-face fluxes separately before taking the conservative two-dimensional divergence.
Field2D compute_rhs(const Field2D& field, const Grid2D& grid, const RuntimeConfig& config) {
    const double hx = dx(grid);
    const double hy = dy(grid);

    std::vector<Conservative2D> fx(static_cast<std::size_t>((grid.nx + 1) * grid.ny));
    std::vector<Conservative2D> fy(static_cast<std::size_t>(grid.nx * (grid.ny + 1)));

    auto fx_at = [&](int i, int j) -> Conservative2D& {
        return fx[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx + 1) + static_cast<std::size_t>(i)];
    };
    auto fy_at = [&](int i, int j) -> Conservative2D& {
        return fy[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx) + static_cast<std::size_t>(i)];
    };

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int iface = 0; iface <= grid.nx; ++iface) {
            const Conservative2D left_state = reconstruct_left_state_x(field, iface - 1, j, config);
            const Conservative2D right_state = reconstruct_right_state_x(field, iface - 1, j, config);
            fx_at(iface, j) = rusanov_flux_x(left_state, right_state, config.gamma);
        }
    }

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int jface = 0; jface <= grid.ny; ++jface) {
        for (int i = 0; i < grid.nx; ++i) {
            const Conservative2D bottom_state = reconstruct_left_state_y(field, i, jface - 1, config);
            const Conservative2D top_state = reconstruct_right_state_y(field, i, jface - 1, config);
            fy_at(i, jface) = rusanov_flux_y(bottom_state, top_state, config.gamma);
        }
    }

    Field2D rhs;
    rhs.nx = grid.nx;
    rhs.ny = grid.ny;
    rhs.cells.resize(field.cells.size());

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            rhs(i, j) = (-1.0 / hx) * (fx_at(i + 1, j) - fx_at(i, j))
                      + (-1.0 / hy) * (fy_at(i, j + 1) - fy_at(i, j));
        }
    }

    return rhs;
}

}  // namespace

// OpenMP changes loop execution only; it uses the same reconstruction, flux, and SSPRK2 update as serial.
SimulationResult2D run_weno2_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state) {
    Field2D field = initial_state;

    double time = 0.0;
    int step = 0;
    while (time < config.final_time && step < config.max_steps) {
        double dt = compute_time_step(field, grid, config);
        dt = std::min(dt, config.final_time - time);

        const Field2D rhs_1 = compute_rhs(field, grid, config);
        const Field2D stage_1 = add_scaled(field, rhs_1, dt);
        const Field2D rhs_2 = compute_rhs(stage_1, grid, config);

        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (std::size_t k = 0; k < field.cells.size(); ++k) {
            field.cells[k] = 0.5 * field.cells[k] + 0.5 * (stage_1.cells[k] + dt * rhs_2.cells[k]);
        }

        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field = std::move(field);
    result.final_time = time;
    result.steps = step;
    return result;
}

}  // namespace fp
