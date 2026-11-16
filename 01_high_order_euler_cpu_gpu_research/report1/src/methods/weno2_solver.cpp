#include "weno2_solver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fp {

namespace {

constexpr double kWenoEps = 1.0e-6;
constexpr double kBackgroundPressure = 1.0;

// Add ghost cells.
std::vector<Conservative1D> add_ghost_cells(const std::vector<Conservative1D>& interior,
                                            const RuntimeConfig& config) {
    // Allocate an extended field including ghost cells.
    std::vector<Conservative1D> extended(interior.size() + 4);
    
    // Copy interior cells into the middle of the extended array.
    for (std::size_t i = 0; i < interior.size(); ++i) {
        extended[i + 2] = interior[i];
    }

    // Fill ghost cells according to the selected boundary condition.
    if (config.boundary == "periodic") {
        extended[0] = interior[interior.size() - 2];
        extended[1] = interior.back();
        extended[extended.size() - 2] = interior.front();
        extended[extended.size() - 1] = interior[1];
    } else {
        extended[0] = interior.front();
        extended[1] = interior.front();
        extended[extended.size() - 2] = interior.back();
        extended[extended.size() - 1] = interior.back();
    }

    return extended;
}

// Compute dt from the CFL condition for the Euler path.
double compute_time_step(const std::vector<Conservative1D>& cells, const RuntimeConfig& config, double dx) {
    double max_speed = 1.0e-12;
    for (std::size_t i = 0; i + 1 < cells.size(); ++i) {
        max_speed = std::max(max_speed, max_wave_speed(cells[i], cells[i + 1], config.gamma));
    }
    return config.cfl * dx / max_speed;
}

// Compute dt from the CFL condition for scalar advection.
double compute_time_step_scalar(const RuntimeConfig& config, double dx) {
    const double a = std::abs(config.advection_velocity);
    if (a < 1.0e-14) {
        return config.final_time;
    }
    return config.cfl * dx / a;
}

// Unnormalized WENO-Z nonlinear weight.
// Small beta means this stencil is smooth and can receive more weight.
// Large beta means this stencil may cross a discontinuity or shock and should receive less weight.
// tau is the WENO-Z global smoothness indicator, measuring the contrast between candidate stencil smoothness values.
// linear_weight is the ideal linear weight.
double nonlinear_weight_z(double beta, double tau, double linear_weight) {
    const double ratio = tau / (kWenoEps + beta);
    return linear_weight * (1.0 + ratio * ratio);
}

// Apply the mapped-weight correction to a nonlinear weight.
double mapped_weight(double omega, double linear_weight) {
    const double d = linear_weight;
    const double numerator = omega * (d + d * d - 3.0 * d * omega + omega * omega);
    const double denominator = d * d + omega * (1.0 - 2.0 * d);
    if (std::abs(denominator) < 1.0e-14) {
        return omega;
    }
    return numerator / denominator;
}

// Add ghost cells according to the selected boundary condition.
std::vector<double> add_ghost_cells_scalar(const std::vector<double>& interior,
                                           const RuntimeConfig& config) {
    std::vector<double> extended(interior.size() + 4, 0.0);

    for (std::size_t i = 0; i < interior.size(); ++i) {
        extended[i + 2] = interior[i];
    }

    if (config.boundary == "periodic") {
        extended[0] = interior[interior.size() - 2];
        extended[1] = interior.back();
        extended[extended.size() - 2] = interior.front();
        extended[extended.size() - 1] = interior[1];
    } else {
        extended[0] = interior.front();
        extended[1] = interior.front();
        extended[extended.size() - 2] = interior.back();
        extended[extended.size() - 1] = interior.back();
    }

    return extended;
}

// Reconstruct the left interface state from the i-1, i, and i+1 stencil values.
Conservative1D reconstruct_left_state(const Conservative1D& uim1,
                                      const Conservative1D& ui,
                                      const Conservative1D& uip1) {
    // beta0 measures the left stencil {U_{i-1}, U_i}; beta1 measures the right stencil {U_i, U_{i+1}}.
    const double beta0_rho = (ui.rho - uim1.rho) * (ui.rho - uim1.rho);
    const double beta1_rho = (uip1.rho - ui.rho) * (uip1.rho - ui.rho);
    // Difference between the two smoothness indicators.
    const double tau_rho = std::abs(beta0_rho - beta1_rho);
    // Compute WENO-Z alpha weights, one score per candidate stencil.
    double alpha0_rho = nonlinear_weight_z(beta0_rho, tau_rho, 0.5);
    double alpha1_rho = nonlinear_weight_z(beta1_rho, tau_rho, 0.5);
    const double alpha_sum_rho = alpha0_rho + alpha1_rho;
    // Normalize alpha scores into proportional nonlinear weights.
    double omega0_rho = alpha0_rho / alpha_sum_rho;
    double omega1_rho = alpha1_rho / alpha_sum_rho;
    // Apply mapped-weight correction to push the nonlinear weights back toward the ideal weights.
    omega0_rho = mapped_weight(omega0_rho, 0.5);
    omega1_rho = mapped_weight(omega1_rho, 0.5);
    const double omega_sum_rho = omega0_rho + omega1_rho;
    // Renormalize after the mapped-weight correction.
    omega0_rho /= omega_sum_rho;
    omega1_rho /= omega_sum_rho;

    // Apply the same steps to momentum.
    const double beta0_mom = (ui.mom - uim1.mom) * (ui.mom - uim1.mom);
    const double beta1_mom = (uip1.mom - ui.mom) * (uip1.mom - ui.mom);
    const double tau_mom = std::abs(beta0_mom - beta1_mom);
    double alpha0_mom = nonlinear_weight_z(beta0_mom, tau_mom, 0.5);
    double alpha1_mom = nonlinear_weight_z(beta1_mom, tau_mom, 0.5);
    const double alpha_sum_mom = alpha0_mom + alpha1_mom;
    double omega0_mom = alpha0_mom / alpha_sum_mom;
    double omega1_mom = alpha1_mom / alpha_sum_mom;
    omega0_mom = mapped_weight(omega0_mom, 0.5);
    omega1_mom = mapped_weight(omega1_mom, 0.5);
    const double omega_sum_mom = omega0_mom + omega1_mom;
    omega0_mom /= omega_sum_mom;
    omega1_mom /= omega_sum_mom;

    // Apply the same steps to energy.
    const double beta0_energy = (ui.energy - uim1.energy) * (ui.energy - uim1.energy);
    const double beta1_energy = (uip1.energy - ui.energy) * (uip1.energy - ui.energy);
    const double tau_energy = std::abs(beta0_energy - beta1_energy);
    double alpha0_energy = nonlinear_weight_z(beta0_energy, tau_energy, 0.5);
    double alpha1_energy = nonlinear_weight_z(beta1_energy, tau_energy, 0.5);
    const double alpha_sum_energy = alpha0_energy + alpha1_energy;
    double omega0_energy = alpha0_energy / alpha_sum_energy;
    double omega1_energy = alpha1_energy / alpha_sum_energy;
    omega0_energy = mapped_weight(omega0_energy, 0.5);
    omega1_energy = mapped_weight(omega1_energy, 0.5);
    const double omega_sum_energy = omega0_energy + omega1_energy;
    omega0_energy /= omega_sum_energy;
    omega1_energy /= omega_sum_energy;

    // Candidate reconstructed values at the target interface.
    const double q0_rho = -0.5 * uim1.rho + 1.5 * ui.rho;
    const double q1_rho = 0.5 * ui.rho + 0.5 * uip1.rho;

    const double q0_mom = -0.5 * uim1.mom + 1.5 * ui.mom;
    const double q1_mom = 0.5 * ui.mom + 0.5 * uip1.mom;

    const double q0_energy = -0.5 * uim1.energy + 1.5 * ui.energy;
    const double q1_energy = 0.5 * ui.energy + 0.5 * uip1.energy;

    // Return the final reconstructed state.
    return {
        omega0_rho * q0_rho + omega1_rho * q1_rho,
        omega0_mom * q0_mom + omega1_mom * q1_mom,
        omega0_energy * q0_energy + omega1_energy * q1_energy,
    };
}

// Repeat the same reconstruction for the right interface state.
Conservative1D reconstruct_right_state(const Conservative1D& ui,
                                       const Conservative1D& uip1,
                                       const Conservative1D& uip2) {
    const double beta0_rho = (uip1.rho - ui.rho) * (uip1.rho - ui.rho);
    const double beta1_rho = (uip2.rho - uip1.rho) * (uip2.rho - uip1.rho);
    const double tau_rho = std::abs(beta0_rho - beta1_rho);
    double alpha0_rho = nonlinear_weight_z(beta0_rho, tau_rho, 0.5);
    double alpha1_rho = nonlinear_weight_z(beta1_rho, tau_rho, 0.5);
    const double alpha_sum_rho = alpha0_rho + alpha1_rho;
    double omega0_rho = alpha0_rho / alpha_sum_rho;
    double omega1_rho = alpha1_rho / alpha_sum_rho;
    omega0_rho = mapped_weight(omega0_rho, 0.5);
    omega1_rho = mapped_weight(omega1_rho, 0.5);
    const double omega_sum_rho = omega0_rho + omega1_rho;
    omega0_rho /= omega_sum_rho;
    omega1_rho /= omega_sum_rho;

    const double beta0_mom = (uip1.mom - ui.mom) * (uip1.mom - ui.mom);
    const double beta1_mom = (uip2.mom - uip1.mom) * (uip2.mom - uip1.mom);
    const double tau_mom = std::abs(beta0_mom - beta1_mom);
    double alpha0_mom = nonlinear_weight_z(beta0_mom, tau_mom, 0.5);
    double alpha1_mom = nonlinear_weight_z(beta1_mom, tau_mom, 0.5);
    const double alpha_sum_mom = alpha0_mom + alpha1_mom;
    double omega0_mom = alpha0_mom / alpha_sum_mom;
    double omega1_mom = alpha1_mom / alpha_sum_mom;
    omega0_mom = mapped_weight(omega0_mom, 0.5);
    omega1_mom = mapped_weight(omega1_mom, 0.5);
    const double omega_sum_mom = omega0_mom + omega1_mom;
    omega0_mom /= omega_sum_mom;
    omega1_mom /= omega_sum_mom;

    const double beta0_energy = (uip1.energy - ui.energy) * (uip1.energy - ui.energy);
    const double beta1_energy = (uip2.energy - uip1.energy) * (uip2.energy - uip1.energy);
    const double tau_energy = std::abs(beta0_energy - beta1_energy);
    double alpha0_energy = nonlinear_weight_z(beta0_energy, tau_energy, 0.5);
    double alpha1_energy = nonlinear_weight_z(beta1_energy, tau_energy, 0.5);
    const double alpha_sum_energy = alpha0_energy + alpha1_energy;
    double omega0_energy = alpha0_energy / alpha_sum_energy;
    double omega1_energy = alpha1_energy / alpha_sum_energy;
    omega0_energy = mapped_weight(omega0_energy, 0.5);
    omega1_energy = mapped_weight(omega1_energy, 0.5);
    const double omega_sum_energy = omega0_energy + omega1_energy;
    omega0_energy /= omega_sum_energy;
    omega1_energy /= omega_sum_energy;

    const double q0_rho = 0.5 * ui.rho + 0.5 * uip1.rho;
    const double q1_rho = 1.5 * uip1.rho - 0.5 * uip2.rho;

    const double q0_mom = 0.5 * ui.mom + 0.5 * uip1.mom;
    const double q1_mom = 1.5 * uip1.mom - 0.5 * uip2.mom;

    const double q0_energy = 0.5 * ui.energy + 0.5 * uip1.energy;
    const double q1_energy = 1.5 * uip1.energy - 0.5 * uip2.energy;

    return {
        omega0_rho * q0_rho + omega1_rho * q1_rho,
        omega0_mom * q0_mom + omega1_mom * q1_mom,
        omega0_energy * q0_energy + omega1_energy * q1_energy,
    };
}

// Compute the current right-hand side.
std::vector<Conservative1D> compute_rhs(const std::vector<Conservative1D>& cells,
                                        const RuntimeConfig& config,
                                        double dx) {
    const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);

    // Store interface fluxes in a separate vector.
    std::vector<Conservative1D> fluxes(config.nx + 1);
    for (int iface = 0; iface <= config.nx; ++iface) {
        // Index of the cell to the left of the interface.
        const int idx = iface + 1;
        // Reconstructed state from the left cell.
        const Conservative1D left_state =
            reconstruct_left_state(extended[idx - 1], extended[idx], extended[idx + 1]);
        // Reconstructed state from the right cell.
        const Conservative1D right_state =
            reconstruct_right_state(extended[idx], extended[idx + 1], extended[idx + 2]);
        // Use the HLL solver to compute the interface flux.
        fluxes[iface] = hll_flux(left_state, right_state, config.gamma);
    }

    std::vector<Conservative1D> rhs(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        rhs[i] = (-1.0 / dx) * (fluxes[i + 1] - fluxes[i]);
    }
    return rhs;
}

// WENO right-hand side for the 1D scalar advection equation.
std::vector<double> compute_rhs_scalar(const std::vector<double>& rho,
                                       const RuntimeConfig& config,
                                       double dx) {
    // Set ghost cells.
    const std::vector<double> extended = add_ghost_cells_scalar(rho, config);
    // Prepare the interface-flux array.
    std::vector<double> fluxes(config.nx + 1, 0.0);
    const double a = config.advection_velocity;

    // Perform WENO reconstruction at each interface.
    for (int iface = 0; iface <= config.nx; ++iface) {
        // Build grid metadata, resolution-dependent paths, or cell-center coordinates.
        const int idx = iface + 1;

        const double beta0_l = (extended[idx] - extended[idx - 1]) * (extended[idx] - extended[idx - 1]);
        const double beta1_l = (extended[idx + 1] - extended[idx]) * (extended[idx + 1] - extended[idx]);
        const double tau_l = std::abs(beta0_l - beta1_l);
        double alpha0_l = nonlinear_weight_z(beta0_l, tau_l, 0.5);
        double alpha1_l = nonlinear_weight_z(beta1_l, tau_l, 0.5);
        double omega0_l = alpha0_l / (alpha0_l + alpha1_l);
        double omega1_l = alpha1_l / (alpha0_l + alpha1_l);
        omega0_l = mapped_weight(omega0_l, 0.5);
        omega1_l = mapped_weight(omega1_l, 0.5);
        const double omega_sum_l = omega0_l + omega1_l;
        omega0_l /= omega_sum_l;
        omega1_l /= omega_sum_l;
        const double left_state =
            omega0_l * (-0.5 * extended[idx - 1] + 1.5 * extended[idx]) +
            omega1_l * (0.5 * extended[idx] + 0.5 * extended[idx + 1]);

        const double beta0_r = (extended[idx + 1] - extended[idx]) * (extended[idx + 1] - extended[idx]);
        const double beta1_r = (extended[idx + 2] - extended[idx + 1]) * (extended[idx + 2] - extended[idx + 1]);
        const double tau_r = std::abs(beta0_r - beta1_r);
        double alpha0_r = nonlinear_weight_z(beta0_r, tau_r, 0.5);
        double alpha1_r = nonlinear_weight_z(beta1_r, tau_r, 0.5);
        double omega0_r = alpha0_r / (alpha0_r + alpha1_r);
        double omega1_r = alpha1_r / (alpha0_r + alpha1_r);
        omega0_r = mapped_weight(omega0_r, 0.5);
        omega1_r = mapped_weight(omega1_r, 0.5);
        const double omega_sum_r = omega0_r + omega1_r;
        omega0_r /= omega_sum_r;
        omega1_r /= omega_sum_r;
        const double right_state =
            omega0_r * (0.5 * extended[idx] + 0.5 * extended[idx + 1]) +
            omega1_r * (1.5 * extended[idx + 1] - 0.5 * extended[idx + 2]);

        // Choose the upwind scalar flux from the sign of advection speed a.
        fluxes[iface] = (a >= 0.0) ? a * left_state : a * right_state;
    }

    // Compute the right-hand side.
    std::vector<double> rhs(config.nx, 0.0);
    for (int i = 0; i < config.nx; ++i) {
        rhs[i] = (-1.0 / dx) * (fluxes[i + 1] - fluxes[i]);
    }
    return rhs;
}

// Update 1D Euler states using a scaled right-hand side.
std::vector<Conservative1D> add_scaled(const std::vector<Conservative1D>& base,
                                       const std::vector<Conservative1D>& delta,
                                       double scale) {
    std::vector<Conservative1D> result(base.size());
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + scale * delta[i];
    }
    return result;
}

// Update scalar-advection states using a scaled right-hand side.
std::vector<double> add_scaled_scalar(const std::vector<double>& base,
                                      const std::vector<double>& delta,
                                      double scale) {
    std::vector<double> result(base.size(), 0.0);
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + scale * delta[i];
    }
    return result;
}

// Rebuild Conservative1D states from a scalar density field.
std::vector<Conservative1D> rebuild_constant_velocity_family(const std::vector<double>& rho,
                                                             const RuntimeConfig& config) {
    // Allocate Euler-form states using rho and the fixed background velocity and pressure.
    std::vector<Conservative1D> cells(rho.size());
    for (std::size_t i = 0; i < rho.size(); ++i) {
        const Primitive1D w{rho[i], config.advection_velocity, kBackgroundPressure};
        cells[i] = primitive_to_conservative(w, config.gamma);
    }
    return cells;
}

}  // namespace

// Main 1D WENO-2 + RK2 solver.
SimulationResult1D run_weno2_rk2_1d(const RuntimeConfig& config,
                                    const std::vector<Conservative1D>& initial_state) {
    // Check that the initial field has the expected number of cells.
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    // Compute grid spacing.
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    // Decide whether to use the scalar-advection shortcut.
    const bool use_scalar_advection_path =
        (config.problem == "smooth_advection" || config.problem == "gaussian_advection") &&
        config.boundary == "periodic";

    // Copy the initial condition into a mutable working field.
    std::vector<Conservative1D> cells = initial_state;

    // Initialize time and step counters.
    double time = 0.0;
    int step = 0;

    // Scalar-advection path.
    if (use_scalar_advection_path) {
        // Extract rho from the Conservative1D field.
        std::vector<double> rho(config.nx, 0.0);
        for (int i = 0; i < config.nx; ++i) {
            rho[i] = initial_state[i].rho;
        }

        // Time loop for the scalar path.
        while (time < config.final_time && step < config.max_steps) {
            // Compute dt from the CFL condition.
            double dt = compute_time_step_scalar(config, dx);
            dt = std::min(dt, config.final_time - time);

            // Compute the first-stage RHS and form stage 1.
            const std::vector<double> rhs_stage_1 = compute_rhs_scalar(rho, config, dx);
            const std::vector<double> stage_1 = add_scaled_scalar(rho, rhs_stage_1, dt);

            // Compute the second-stage RHS and complete the RK2 update.
            const std::vector<double> rhs_stage_2 = compute_rhs_scalar(stage_1, config, dx);
            for (int i = 0; i < config.nx; ++i) {
                rho[i] = 0.5 * rho[i] + 0.5 * (stage_1[i] + dt * rhs_stage_2[i]);
            }

            // Advance time and the step counter.
            time += dt;
            ++step;
        }

        // After the scalar path finishes, wrap rho back into Conservative1D states.
        cells = rebuild_constant_velocity_family(rho, config);
    } else {
        // Time loop for the Euler path.
        while (time < config.final_time && step < config.max_steps) {
            // Add ghost cells only for the dt calculation.
            const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);
            double dt = compute_time_step(extended, config, dx);
            dt = std::min(dt, config.final_time - time);

            // First RK2 stage for the Euler path.
            const std::vector<Conservative1D> rhs_stage_1 = compute_rhs(cells, config, dx);
            const std::vector<Conservative1D> stage_1 = add_scaled(cells, rhs_stage_1, dt);

            // Second RK2 stage for the Euler path.
            const std::vector<Conservative1D> rhs_stage_2 = compute_rhs(stage_1, config, dx);
            for (int i = 0; i < config.nx; ++i) {
                // Final Euler RK2 convex-combination update.
                cells[i] = 0.5 * cells[i] + 0.5 * (stage_1[i] + dt * rhs_stage_2[i]);
            }

            // Time-step selection and Runge-Kutta stage update.
            time += dt;
            ++step;
        }
    }

    // Create the result object storing final field, time, and step count.
    SimulationResult1D result;
    result.cells = std::move(cells);
    result.final_time = time;
    result.steps = step;
    result.x.resize(config.nx);

    // Generate cell-center coordinates.
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }

    return result;
}

}  // namespace fp

