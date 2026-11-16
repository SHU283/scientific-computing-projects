#include "weno3_solver_2d.h"

#include "phase_profiler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
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

int wrap_index(int index, int size) {
    index %= size;
    return index < 0 ? index + size : index;
}

const Conservative2D& sample_cell(const Field2D& field,
                                  int i,
                                  int j,
                                  const RuntimeConfig& config) {
    if (config.boundary == "periodic") {
        return field(wrap_index(i, field.nx), wrap_index(j, field.ny));
    }
    return field(std::clamp(i, 0, field.nx - 1),
                 std::clamp(j, 0, field.ny - 1));
}

double nonlinear_weight_z(double beta, double tau, double linear_weight) {
    const double ratio = tau / (kWenoEps + beta);
    return linear_weight * (1.0 + ratio * ratio);
}

double mapped_weight(double omega, double linear_weight) {
    const double d = linear_weight;
    const double numerator = omega * (d + d * d - 3.0 * d * omega + omega * omega);
    const double denominator = d * d + omega * (1.0 - 2.0 * d);
    return std::abs(denominator) < 1.0e-14 ? omega : numerator / denominator;
}

// Reconstruct each conservative component before forming directional HLLC fluxes.
double reconstruct_left_component(double uim1, double ui, double uip1) {
    const double q0 = -0.5 * uim1 + 1.5 * ui;
    const double q1 = 0.5 * ui + 0.5 * uip1;
    const double beta0 = (ui - uim1) * (ui - uim1);
    const double beta1 = (uip1 - ui) * (uip1 - ui);
    const double tau = std::abs(beta0 - beta1);

    double omega0 = nonlinear_weight_z(beta0, tau, 1.0 / 3.0);
    double omega1 = nonlinear_weight_z(beta1, tau, 2.0 / 3.0);
    const double alpha_sum = omega0 + omega1;
    omega0 = mapped_weight(omega0 / alpha_sum, 1.0 / 3.0);
    omega1 = mapped_weight(omega1 / alpha_sum, 2.0 / 3.0);
    return (omega0 * q0 + omega1 * q1) / (omega0 + omega1);
}

double reconstruct_right_component(double ui, double uip1, double uip2) {
    const double q0 = 0.5 * ui + 0.5 * uip1;
    const double q1 = 1.5 * uip1 - 0.5 * uip2;
    const double beta0 = (uip1 - ui) * (uip1 - ui);
    const double beta1 = (uip2 - uip1) * (uip2 - uip1);
    const double tau = std::abs(beta0 - beta1);

    double omega0 = nonlinear_weight_z(beta0, tau, 2.0 / 3.0);
    double omega1 = nonlinear_weight_z(beta1, tau, 1.0 / 3.0);
    const double alpha_sum = omega0 + omega1;
    omega0 = mapped_weight(omega0 / alpha_sum, 2.0 / 3.0);
    omega1 = mapped_weight(omega1 / alpha_sum, 1.0 / 3.0);
    return (omega0 * q0 + omega1 * q1) / (omega0 + omega1);
}

Conservative2D reconstruct_left_x(const Field2D& field,
                                  int i,
                                  int j,
                                  const RuntimeConfig& config) {
    const Conservative2D& uim1 = sample_cell(field, i - 1, j, config);
    const Conservative2D& ui = sample_cell(field, i, j, config);
    const Conservative2D& uip1 = sample_cell(field, i + 1, j, config);
    return {
        reconstruct_left_component(uim1.rho, ui.rho, uip1.rho),
        reconstruct_left_component(uim1.mom_x, ui.mom_x, uip1.mom_x),
        reconstruct_left_component(uim1.mom_y, ui.mom_y, uip1.mom_y),
        reconstruct_left_component(uim1.energy, ui.energy, uip1.energy),
    };
}

Conservative2D reconstruct_right_x(const Field2D& field,
                                   int i,
                                   int j,
                                   const RuntimeConfig& config) {
    const Conservative2D& ui = sample_cell(field, i, j, config);
    const Conservative2D& uip1 = sample_cell(field, i + 1, j, config);
    const Conservative2D& uip2 = sample_cell(field, i + 2, j, config);
    return {
        reconstruct_right_component(ui.rho, uip1.rho, uip2.rho),
        reconstruct_right_component(ui.mom_x, uip1.mom_x, uip2.mom_x),
        reconstruct_right_component(ui.mom_y, uip1.mom_y, uip2.mom_y),
        reconstruct_right_component(ui.energy, uip1.energy, uip2.energy),
    };
}

Conservative2D reconstruct_left_y(const Field2D& field,
                                  int i,
                                  int j,
                                  const RuntimeConfig& config) {
    const Conservative2D& ujm1 = sample_cell(field, i, j - 1, config);
    const Conservative2D& uj = sample_cell(field, i, j, config);
    const Conservative2D& ujp1 = sample_cell(field, i, j + 1, config);
    return {
        reconstruct_left_component(ujm1.rho, uj.rho, ujp1.rho),
        reconstruct_left_component(ujm1.mom_x, uj.mom_x, ujp1.mom_x),
        reconstruct_left_component(ujm1.mom_y, uj.mom_y, ujp1.mom_y),
        reconstruct_left_component(ujm1.energy, uj.energy, ujp1.energy),
    };
}

Conservative2D reconstruct_right_y(const Field2D& field,
                                   int i,
                                   int j,
                                   const RuntimeConfig& config) {
    const Conservative2D& uj = sample_cell(field, i, j, config);
    const Conservative2D& ujp1 = sample_cell(field, i, j + 1, config);
    const Conservative2D& ujp2 = sample_cell(field, i, j + 2, config);
    return {
        reconstruct_right_component(uj.rho, ujp1.rho, ujp2.rho),
        reconstruct_right_component(uj.mom_x, ujp1.mom_x, ujp2.mom_x),
        reconstruct_right_component(uj.mom_y, ujp1.mom_y, ujp2.mom_y),
        reconstruct_right_component(uj.energy, ujp1.energy, ujp2.energy),
    };
}

bool is_admissible(const Conservative2D& state, double gamma) {
    if (!(state.rho > 0.0) || !std::isfinite(state.rho) ||
        !std::isfinite(state.mom_x) || !std::isfinite(state.mom_y) ||
        !std::isfinite(state.energy)) {
        return false;
    }
    const double kinetic =
        0.5 * (state.mom_x * state.mom_x + state.mom_y * state.mom_y) / state.rho;
    return (gamma - 1.0) * (state.energy - kinetic) > 0.0;
}

double compute_time_step(const Field2D& field,
                         const Grid2D& grid,
                         const RuntimeConfig& config) {
    double max_x = 1.0e-12;
    double max_y = 1.0e-12;
    #pragma omp parallel for collapse(2) reduction(max:max_x,max_y) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < field.ny; ++j) {
        for (int i = 0; i < field.nx; ++i) {
            max_x = std::max(max_x, max_signal_speed_x(field(i, j), config.gamma));
            max_y = std::max(max_y, max_signal_speed_y(field(i, j), config.gamma));
        }
    }
    return config.cfl / (max_x / dx(grid) + max_y / dy(grid));
}

Field2D add_scaled(const Field2D& base,
                   const Field2D& delta,
                   double scale,
                   const RuntimeConfig& config) {
    Field2D result = base;
    #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
    for (std::size_t k = 0; k < base.cells.size(); ++k) {
        result.cells[k] = base.cells[k] + scale * delta.cells[k];
    }
    return result;
}

Field2D compute_rhs(const Field2D& field,
                    const Grid2D& grid,
                    const RuntimeConfig& config,
                    PhaseProfile* profile) {
    const auto workspace_start = PhaseProfile::Clock::now();
    std::vector<Conservative2D> flux_x(
        static_cast<std::size_t>(grid.nx + 1) * static_cast<std::size_t>(grid.ny));
    std::vector<Conservative2D> flux_y(
        static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny + 1));

    auto x_flux = [&](int i, int j) -> Conservative2D& {
        return flux_x[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx + 1) +
                      static_cast<std::size_t>(i)];
    };
    auto y_flux = [&](int i, int j) -> Conservative2D& {
        return flux_y[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx) +
                      static_cast<std::size_t>(i)];
    };
    if (profile) {
        profile->workspace_s += PhaseProfile::elapsed(workspace_start);
    }

    const auto flux_x_start = PhaseProfile::Clock::now();
    // OpenMP changes loop scheduling only; reconstruction and flux formulas are shared.
    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int iface = 0; iface <= grid.nx; ++iface) {
            Conservative2D left = reconstruct_left_x(field, iface - 1, j, config);
            Conservative2D right = reconstruct_right_x(field, iface - 1, j, config);
            if (!is_admissible(left, config.gamma)) {
                left = sample_cell(field, iface - 1, j, config);
            }
            if (!is_admissible(right, config.gamma)) {
                right = sample_cell(field, iface, j, config);
            }
            x_flux(iface, j) = hllc_flux_x(left, right, config.gamma);
        }
    }
    if (profile) {
        profile->flux_x_s += PhaseProfile::elapsed(flux_x_start);
    }

    const auto flux_y_start = PhaseProfile::Clock::now();
    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int jface = 0; jface <= grid.ny; ++jface) {
        for (int i = 0; i < grid.nx; ++i) {
            Conservative2D bottom = reconstruct_left_y(field, i, jface - 1, config);
            Conservative2D top = reconstruct_right_y(field, i, jface - 1, config);
            if (!is_admissible(bottom, config.gamma)) {
                bottom = sample_cell(field, i, jface - 1, config);
            }
            if (!is_admissible(top, config.gamma)) {
                top = sample_cell(field, i, jface, config);
            }
            y_flux(i, jface) = hllc_flux_y(bottom, top, config.gamma);
        }
    }
    if (profile) {
        profile->flux_y_s += PhaseProfile::elapsed(flux_y_start);
    }

    const auto rhs_start = PhaseProfile::Clock::now();
    Field2D rhs;
    rhs.nx = grid.nx;
    rhs.ny = grid.ny;
    rhs.cells.resize(field.cells.size());
    const double inv_dx = 1.0 / dx(grid);
    const double inv_dy = 1.0 / dy(grid);

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            rhs(i, j) = -inv_dx * (x_flux(i + 1, j) - x_flux(i, j))
                      - inv_dy * (y_flux(i, j + 1) - y_flux(i, j));
        }
    }
    if (profile) {
        profile->rhs_assembly_s += PhaseProfile::elapsed(rhs_start);
    }
    return rhs;
}

}  // namespace

SimulationResult2D run_weno3_rk3_2d(const RuntimeConfig& config,
                                    const Grid2D& grid,
                                    const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("Initial 2D field does not match the configured grid.");
    }

    PhaseProfile profile;
    profile.enabled = PhaseProfile::requested();
    const auto total_start = PhaseProfile::Clock::now();
    const auto setup_start = PhaseProfile::Clock::now();
    Field2D field = initial_state;
    if (profile.enabled) {
        profile.setup_s += PhaseProfile::elapsed(setup_start);
    }
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        const auto cfl_start = PhaseProfile::Clock::now();
        double dt = compute_time_step(field, grid, config);
        dt = std::min(dt, config.final_time - time);
        if (profile.enabled) {
            profile.cfl_s += PhaseProfile::elapsed(cfl_start);
        }

        // SSPRK3 combines three evaluations of the same conservative residual.
        const Field2D rhs1 = compute_rhs(
            field, grid, config, profile.enabled ? &profile : nullptr);
        const auto update_1_start = PhaseProfile::Clock::now();
        const Field2D stage1 = add_scaled(field, rhs1, dt, config);
        if (profile.enabled) {
            profile.state_update_s += PhaseProfile::elapsed(update_1_start);
        }
        const Field2D rhs2 = compute_rhs(
            stage1, grid, config, profile.enabled ? &profile : nullptr);

        const auto update_2_start = PhaseProfile::Clock::now();
        Field2D stage2 = field;
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (std::size_t k = 0; k < field.cells.size(); ++k) {
            stage2.cells[k] =
                0.75 * field.cells[k] + 0.25 * (stage1.cells[k] + dt * rhs2.cells[k]);
        }
        if (profile.enabled) {
            profile.state_update_s += PhaseProfile::elapsed(update_2_start);
        }

        const Field2D rhs3 = compute_rhs(
            stage2, grid, config, profile.enabled ? &profile : nullptr);
        const auto update_3_start = PhaseProfile::Clock::now();
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (std::size_t k = 0; k < field.cells.size(); ++k) {
            field.cells[k] =
                (1.0 / 3.0) * field.cells[k] +
                (2.0 / 3.0) * (stage2.cells[k] + dt * rhs3.cells[k]);
        }
        if (profile.enabled) {
            profile.state_update_s += PhaseProfile::elapsed(update_3_start);
        }

        time += dt;
        ++step;
    }

    profile.total_s = PhaseProfile::elapsed(total_start);
    profile.print("weno3", config.backend, grid.nx, grid.ny, step);
    return {grid, std::move(field), time, step};
}

}  // namespace fp
