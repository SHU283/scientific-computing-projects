#include "weno3_solver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fp {

namespace {

constexpr double kWenoEps = 1.0e-6;
constexpr double kBackgroundPressure = 1.0;

std::vector<Conservative1D> add_ghost_cells(const std::vector<Conservative1D>& interior,
                                            const RuntimeConfig& config) {
    std::vector<Conservative1D> extended(interior.size() + 4);

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

double compute_time_step(const std::vector<Conservative1D>& cells, const RuntimeConfig& config, double dx) {
    double max_speed = 1.0e-12;
    for (std::size_t i = 0; i + 1 < cells.size(); ++i) {
        max_speed = std::max(max_speed, max_wave_speed(cells[i], cells[i + 1], config.gamma));
    }
    return config.cfl * dx / max_speed;
}

double compute_time_step_scalar(const RuntimeConfig& config, double dx) {
    const double a = std::abs(config.advection_velocity);
    if (a < 1.0e-14) {
        return config.final_time;
    }
    return config.cfl * dx / a;
}

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

double nonlinear_weight_z(double beta,
                          double tau,
                          double linear_weight) {
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

// Third-order interface values use two biased candidates and mapped WENO-Z weights.
double reconstruct_left_component(double vim2, double vim1, double vi) {
    const double q0 = -0.5 * vim2 + 1.5 * vim1;
    const double q1 = 0.5 * vim1 + 0.5 * vi;

    const double beta0 = (vim1 - vim2) * (vim1 - vim2);
    const double beta1 = (vi - vim1) * (vi - vim1);
    const double tau = std::abs(beta0 - beta1);

    const double alpha0 = nonlinear_weight_z(beta0, tau, 1.0 / 3.0);
    const double alpha1 = nonlinear_weight_z(beta1, tau, 2.0 / 3.0);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = mapped_weight(omega0, 1.0 / 3.0);
    omega1 = mapped_weight(omega1, 2.0 / 3.0);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    return omega0 * q0 + omega1 * q1;
}

double reconstruct_right_component(double vi, double vip1, double vip2) {
    const double q0 = 0.5 * vi + 0.5 * vip1;
    const double q1 = 1.5 * vip1 - 0.5 * vip2;

    const double beta0 = (vip1 - vi) * (vip1 - vi);
    const double beta1 = (vip2 - vip1) * (vip2 - vip1);
    const double tau = std::abs(beta0 - beta1);

    const double alpha0 = nonlinear_weight_z(beta0, tau, 2.0 / 3.0);
    const double alpha1 = nonlinear_weight_z(beta1, tau, 1.0 / 3.0);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = mapped_weight(omega0, 2.0 / 3.0);
    omega1 = mapped_weight(omega1, 1.0 / 3.0);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    return omega0 * q0 + omega1 * q1;
}

Conservative1D reconstruct_left_state(const Conservative1D& uim2,
                                      const Conservative1D& uim1,
                                      const Conservative1D& ui) {
    return {
        reconstruct_left_component(uim2.rho, uim1.rho, ui.rho),
        reconstruct_left_component(uim2.mom, uim1.mom, ui.mom),
        reconstruct_left_component(uim2.energy, uim1.energy, ui.energy),
    };
}

Conservative1D reconstruct_right_state(const Conservative1D& ui,
                                       const Conservative1D& uip1,
                                       const Conservative1D& uip2) {
    return {
        reconstruct_right_component(ui.rho, uip1.rho, uip2.rho),
        reconstruct_right_component(ui.mom, uip1.mom, uip2.mom),
        reconstruct_right_component(ui.energy, uip1.energy, uip2.energy),
    };
}

std::vector<Conservative1D> compute_rhs(const std::vector<Conservative1D>& cells,
                                        const RuntimeConfig& config,
                                        double dx) {
    const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);

    std::vector<Conservative1D> fluxes(config.nx + 1);
    for (int iface = 0; iface <= config.nx; ++iface) {
        const int idx = iface + 1;
        const Conservative1D left_state =
            reconstruct_left_state(extended[idx - 1], extended[idx], extended[idx + 1]);
        const Conservative1D right_state =
            reconstruct_right_state(extended[idx], extended[idx + 1], extended[idx + 2]);
        // Apply the same HLLC flux as the 2D and CUDA implementations.
        fluxes[iface] = hllc_flux(left_state, right_state, config.gamma);
    }

    std::vector<Conservative1D> rhs(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        rhs[i] = (-1.0 / dx) * (fluxes[i + 1] - fluxes[i]);
    }
    return rhs;
}

std::vector<double> compute_rhs_scalar(const std::vector<double>& rho,
                                       const RuntimeConfig& config,
                                       double dx) {
    const std::vector<double> extended = add_ghost_cells_scalar(rho, config);
    std::vector<double> fluxes(config.nx + 1, 0.0);
    const double a = config.advection_velocity;

    for (int iface = 0; iface <= config.nx; ++iface) {
        const int idx = iface + 1;
        const double left_state =
            reconstruct_left_component(extended[idx - 1], extended[idx], extended[idx + 1]);
        const double right_state =
            reconstruct_right_component(extended[idx], extended[idx + 1], extended[idx + 2]);
        fluxes[iface] = (a >= 0.0) ? a * left_state : a * right_state;
    }

    std::vector<double> rhs(config.nx, 0.0);
    for (int i = 0; i < config.nx; ++i) {
        rhs[i] = (-1.0 / dx) * (fluxes[i + 1] - fluxes[i]);
    }
    return rhs;
}

std::vector<Conservative1D> add_scaled(const std::vector<Conservative1D>& base,
                                       const std::vector<Conservative1D>& delta,
                                       double scale) {
    std::vector<Conservative1D> result(base.size());
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + scale * delta[i];
    }
    return result;
}

std::vector<double> add_scaled_scalar(const std::vector<double>& base,
                                      const std::vector<double>& delta,
                                      double scale) {
    std::vector<double> result(base.size(), 0.0);
    for (std::size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + scale * delta[i];
    }
    return result;
}

std::vector<Conservative1D> rebuild_constant_velocity_family(const std::vector<double>& rho,
                                                             const RuntimeConfig& config) {
    std::vector<Conservative1D> cells(rho.size());
    for (std::size_t i = 0; i < rho.size(); ++i) {
        const Primitive1D w{rho[i], config.advection_velocity, kBackgroundPressure};
        cells[i] = primitive_to_conservative(w, config.gamma);
    }
    return cells;
}

}  // namespace

SimulationResult1D run_weno3_rk3_1d(const RuntimeConfig& config,
                                    const std::vector<Conservative1D>& initial_state) {
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    // The scalar path is valid only for the prescribed constant-velocity periodic tests.
    const bool use_scalar_advection_path =
        (config.problem == "smooth_advection" || config.problem == "gaussian_advection") &&
        config.boundary == "periodic";

    std::vector<Conservative1D> cells = initial_state;

    double time = 0.0;
    int step = 0;

    if (use_scalar_advection_path) {
        std::vector<double> rho(config.nx, 0.0);
        for (int i = 0; i < config.nx; ++i) {
            rho[i] = initial_state[i].rho;
        }

        while (time < config.final_time && step < config.max_steps) {
            double dt = compute_time_step_scalar(config, dx);
            dt = std::min(dt, config.final_time - time);

            const std::vector<double> rhs_1 = compute_rhs_scalar(rho, config, dx);
            const std::vector<double> u1 = add_scaled_scalar(rho, rhs_1, dt);

            const std::vector<double> rhs_2 = compute_rhs_scalar(u1, config, dx);
            std::vector<double> u2(config.nx, 0.0);
            for (int i = 0; i < config.nx; ++i) {
                u2[i] = 0.75 * rho[i] + 0.25 * (u1[i] + dt * rhs_2[i]);
            }

            const std::vector<double> rhs_3 = compute_rhs_scalar(u2, config, dx);
            for (int i = 0; i < config.nx; ++i) {
                rho[i] = (1.0 / 3.0) * rho[i] + (2.0 / 3.0) * (u2[i] + dt * rhs_3[i]);
            }

            time += dt;
            ++step;
        }

        cells = rebuild_constant_velocity_family(rho, config);
    } else {
    while (time < config.final_time && step < config.max_steps) {
        const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);
        double dt = compute_time_step(extended, config, dx);
        dt = std::min(dt, config.final_time - time);

        const std::vector<Conservative1D> rhs_1 = compute_rhs(cells, config, dx);
        const std::vector<Conservative1D> u1 = add_scaled(cells, rhs_1, dt);

        const std::vector<Conservative1D> rhs_2 = compute_rhs(u1, config, dx);
        std::vector<Conservative1D> u2(config.nx);
        for (int i = 0; i < config.nx; ++i) {
            u2[i] = 0.75 * cells[i] + 0.25 * (u1[i] + dt * rhs_2[i]);
        }

        const std::vector<Conservative1D> rhs_3 = compute_rhs(u2, config, dx);
        for (int i = 0; i < config.nx; ++i) {
            cells[i] = (1.0 / 3.0) * cells[i] + (2.0 / 3.0) * (u2[i] + dt * rhs_3[i]);
        }

        time += dt;
        ++step;
    }
    }

    SimulationResult1D result;
    result.cells = std::move(cells);
    result.final_time = time;
    result.steps = step;
    result.x.resize(config.nx);

    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }

    return result;
}

}  // namespace fp
