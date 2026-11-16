#include "ader3_solver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fp {

namespace {

constexpr double kBackgroundPressure = 1.0;

double periodic_coordinate(double xi, int nx) {
    const double period = static_cast<double>(nx);
    while (xi < 0.0) {
        xi += period;
    }
    while (xi >= period) {
        xi -= period;
    }
    return xi;
}

double rho_from_state(const Conservative1D& cell) {
    return cell.rho;
}

std::vector<double> point_values_to_cell_averages(const std::vector<double>& values) {
    const int nx = static_cast<int>(values.size());
    std::vector<double> averages(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        const int im1 = (i - 1 + nx) % nx;
        const int ip1 = (i + 1) % nx;
        averages[i] = values[i] + (values[im1] - 2.0 * values[i] + values[ip1]) / 24.0;
    }
    return averages;
}

double compute_time_step_scalar(const RuntimeConfig& config, double dx) {
    const double a = std::abs(config.advection_velocity);
    if (a < 1.0e-14) {
        return config.final_time;
    }
    return config.cfl * dx / a;
}

double quadratic_reconstruction(double uim1, double ui, double uip1, double xi) {
    const double delta2 = uip1 - 2.0 * ui + uim1;
    const double slope = 0.5 * (uip1 - uim1);
    return ui - delta2 / 24.0 + slope * xi + 0.5 * delta2 * xi * xi;
}

double integrate_quadratic_segment(double uim1,
                                   double ui,
                                   double uip1,
                                   double xi_left,
                                   double xi_right) {
    const double delta2 = uip1 - 2.0 * ui + uim1;
    const double slope = 0.5 * (uip1 - uim1);
    const double constant = ui - delta2 / 24.0;

    const auto primitive = [&](double xi) {
        return constant * xi + 0.5 * slope * xi * xi + (delta2 / 6.0) * xi * xi * xi;
    };

    return primitive(xi_right) - primitive(xi_left);
}

std::vector<double> density_one_step_ader3(const std::vector<double>& rho,
                                           const RuntimeConfig& config,
                                           double dt,
                                           double dx) {
    const int nx = static_cast<int>(rho.size());
    const double a = config.advection_velocity;
    const double nu = std::abs(a) * dt / dx;

    if (nu > 1.0 + 1.0e-12) {
        throw std::runtime_error("ADER-3 prototype requires CFL <= 1 for the scalar predictor.");
    }

    std::vector<double> next(nx, 0.0);

    if (std::abs(a) < 1.0e-14) {
        return rho;
    }

    for (int i = 0; i < nx; ++i) {
        const int im1 = (i - 1 + nx) % nx;
        const int ip1 = (i + 1) % nx;
        const int ip2 = (i + 2) % nx;

        if (a > 0.0) {
            const double flux_right =
                integrate_quadratic_segment(rho[im1], rho[i], rho[ip1], 0.5 - nu, 0.5);
            const double flux_left =
                integrate_quadratic_segment(rho[(im1 - 1 + nx) % nx], rho[im1], rho[i], 0.5 - nu, 0.5);
            next[i] = rho[i] - (flux_right - flux_left);
        } else {
            const double flux_right =
                integrate_quadratic_segment(rho[i], rho[ip1], rho[ip2], -0.5, -0.5 + nu);
            const double flux_left =
                integrate_quadratic_segment(rho[im1], rho[i], rho[ip1], -0.5, -0.5 + nu);
            next[i] = rho[i] + (flux_right - flux_left);
        }
    }

    return next;
}

std::vector<Conservative1D> rebuild_euler_family(const std::vector<double>& rho,
                                                 const RuntimeConfig& config) {
    std::vector<Conservative1D> cells(rho.size());
    for (std::size_t i = 0; i < rho.size(); ++i) {
        const Primitive1D w{rho[i], config.advection_velocity, kBackgroundPressure};
        cells[i] = primitive_to_conservative(w, config.gamma);
    }
    return cells;
}

}  // namespace

SimulationResult1D run_ader3_1d(const RuntimeConfig& config,
                                const std::vector<Conservative1D>& initial_state) {
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    if (config.problem != "smooth_advection" && config.problem != "gaussian_advection") {
        throw std::runtime_error(
            "Current ADER-3 prototype only supports problem=smooth_advection or problem=gaussian_advection.");
    }

    if (config.boundary != "periodic") {
        throw std::runtime_error("Current ADER-3 prototype requires periodic boundary conditions.");
    }

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    std::vector<Conservative1D> cells = initial_state;
    {
        std::vector<double> rho(config.nx);
        for (int i = 0; i < config.nx; ++i) {
            rho[i] = rho_from_state(cells[i]);
        }
        cells = rebuild_euler_family(point_values_to_cell_averages(rho), config);
    }

    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        double dt = compute_time_step_scalar(config, dx);
        dt = std::min(dt, config.final_time - time);

        std::vector<double> rho(config.nx);
        for (int i = 0; i < config.nx; ++i) {
            rho[i] = rho_from_state(cells[i]);
        }

        const std::vector<double> next_rho = density_one_step_ader3(rho, config, dt, dx);
        cells = rebuild_euler_family(next_rho, config);

        time += dt;
        ++step;
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

