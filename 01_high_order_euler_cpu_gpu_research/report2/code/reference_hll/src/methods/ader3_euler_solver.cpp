#include "ader3_euler_solver.h"

#include <algorithm>
#include <array>
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

constexpr double kReferenceDerivativeStep = 1.0e-3;
constexpr double kDirectionalDifferenceScale = 1.0e-6;

int sample_index(int index, int size, const RuntimeConfig& config) {
    if (config.boundary == "periodic") {
        index %= size;
        return index < 0 ? index + size : index;
    }
    return std::clamp(index, 0, size - 1);
}

Conservative1D sample_cell(const std::vector<Conservative1D>& cells,
                           int index,
                           const RuntimeConfig& config) {
    return cells[sample_index(index, static_cast<int>(cells.size()), config)];
}

bool is_admissible(const Conservative1D& state, double gamma) {
    if (!(state.rho > 0.0) || !std::isfinite(state.rho) ||
        !std::isfinite(state.mom) || !std::isfinite(state.energy)) {
        return false;
    }
    const double kinetic = 0.5 * state.mom * state.mom / state.rho;
    return (gamma - 1.0) * (state.energy - kinetic) > 0.0;
}

Conservative1D evaluate_quadratic(const Conservative1D& um1,
                                  const Conservative1D& u0,
                                  const Conservative1D& up1,
                                  double xi) {
    const Conservative1D second_difference = up1 - 2.0 * u0 + um1;
    const Conservative1D slope = 0.5 * (up1 - um1);
    return u0 - (1.0 / 24.0) * second_difference +
           xi * slope + (0.5 * xi * xi) * second_difference;
}

struct LocalPolynomial {
    Conservative1D um1;
    Conservative1D u0;
    Conservative1D up1;

    Conservative1D evaluate(double xi) const {
        return evaluate_quadratic(um1, u0, up1, xi);
    }
};

bool polynomial_is_admissible(const LocalPolynomial& polynomial, double gamma) {
    constexpr std::array<double, 7> check_points = {
        -0.5 - 2.0 * kReferenceDerivativeStep,
        -0.5,
        -0.25,
        0.0,
        0.25,
        0.5,
        0.5 + 2.0 * kReferenceDerivativeStep,
    };
    for (const double xi : check_points) {
        if (!is_admissible(polynomial.evaluate(xi), gamma)) {
            return false;
        }
    }
    return true;
}

bool interface_is_troubled(const Conservative1D& left,
                           const Conservative1D& right,
                           double gamma) {
    const Primitive1D wl = conservative_to_primitive(left, gamma);
    const Primitive1D wr = conservative_to_primitive(right, gamma);
    const double rho_jump =
        std::abs(wr.rho - wl.rho) / std::max({wl.rho, wr.rho, 1.0e-12});
    const double pressure_jump =
        std::abs(wr.p - wl.p) / std::max({wl.p, wr.p, 1.0e-12});
    const double sound_scale =
        std::max(sound_speed(wl, gamma), sound_speed(wr, gamma));
    const double velocity_jump =
        std::abs(wr.u - wl.u) / std::max(sound_scale, 1.0e-12);
    return std::max({rho_jump, pressure_jump, velocity_jump}) > 0.15;
}

LocalPolynomial reconstruct_cell(const std::vector<Conservative1D>& cells,
                                 int i,
                                 const RuntimeConfig& config) {
    return {
        sample_cell(cells, i - 1, config),
        sample_cell(cells, i, config),
        sample_cell(cells, i + 1, config),
    };
}

Conservative1D spatial_time_derivative(const LocalPolynomial& polynomial,
                                       double xi,
                                       double dx,
                                       double gamma) {
    const double h = kReferenceDerivativeStep;
    const Conservative1D flux_plus = physical_flux(polynomial.evaluate(xi + h), gamma);
    const Conservative1D flux_minus = physical_flux(polynomial.evaluate(xi - h), gamma);
    return (-1.0 / (2.0 * h * dx)) * (flux_plus - flux_minus);
}

double max_component_scale(const Conservative1D& state) {
    return std::max({1.0, std::abs(state.rho), std::abs(state.mom), std::abs(state.energy)});
}

Conservative1D flux_jacobian_vector_product(const Conservative1D& state,
                                            const Conservative1D& direction,
                                            double gamma) {
    const double state_scale = max_component_scale(state);
    const double direction_scale = max_component_scale(direction);
    double epsilon = kDirectionalDifferenceScale * state_scale / direction_scale;

    Conservative1D plus = state + epsilon * direction;
    Conservative1D minus = state - epsilon * direction;
    for (int attempt = 0;
         attempt < 8 &&
         (!is_admissible(plus, gamma) || !is_admissible(minus, gamma));
         ++attempt) {
        epsilon *= 0.25;
        plus = state + epsilon * direction;
        minus = state - epsilon * direction;
    }

    if (!is_admissible(plus, gamma) || !is_admissible(minus, gamma)) {
        throw std::runtime_error(
            "ADER-3 local predictor could not form an admissible Jacobian perturbation.");
    }
    return (physical_flux(plus, gamma) - physical_flux(minus, gamma)) /
           (2.0 * epsilon);
}

Conservative1D second_time_derivative(const LocalPolynomial& polynomial,
                                      double xi,
                                      double dx,
                                      double gamma) {
    const double h = kReferenceDerivativeStep;
    const Conservative1D state_plus = polynomial.evaluate(xi + h);
    const Conservative1D state_minus = polynomial.evaluate(xi - h);
    const Conservative1D ut_plus =
        spatial_time_derivative(polynomial, xi + h, dx, gamma);
    const Conservative1D ut_minus =
        spatial_time_derivative(polynomial, xi - h, dx, gamma);

    const Conservative1D a_ut_plus =
        flux_jacobian_vector_product(state_plus, ut_plus, gamma);
    const Conservative1D a_ut_minus =
        flux_jacobian_vector_product(state_minus, ut_minus, gamma);
    return (-1.0 / (2.0 * h * dx)) * (a_ut_plus - a_ut_minus);
}

Conservative1D predict_local_state(const LocalPolynomial& polynomial,
                                   double xi,
                                   double tau,
                                   double dx,
                                   double gamma) {
    const Conservative1D state = polynomial.evaluate(xi);
    const Conservative1D ut = spatial_time_derivative(polynomial, xi, dx, gamma);
    const Conservative1D utt = second_time_derivative(polynomial, xi, dx, gamma);
    return state + tau * ut + (0.5 * tau * tau) * utt;
}

double compute_time_step(const std::vector<Conservative1D>& cells,
                         const RuntimeConfig& config,
                         double dx) {
    double max_speed = 1.0e-12;
    #pragma omp parallel for reduction(max:max_speed) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int i = 0; i <= config.nx; ++i) {
        const Conservative1D left = sample_cell(cells, i - 1, config);
        const Conservative1D right = sample_cell(cells, i, config);
        max_speed = std::max(max_speed, max_wave_speed(left, right, config.gamma));
    }
    return config.cfl * dx / max_speed;
}

Conservative1D time_averaged_interface_flux(
    const std::vector<Conservative1D>& cells,
    int iface,
    const RuntimeConfig& config,
    double dt,
    double dx) {
    constexpr double offset = 0.28867513459481288225;
    const std::array<double, 2> time_nodes = {
        dt * (0.5 - offset),
        dt * (0.5 + offset),
    };

    const int left_cell = iface - 1;
    const int right_cell = iface;
    const LocalPolynomial left_polynomial =
        reconstruct_cell(cells, left_cell, config);
    const LocalPolynomial right_polynomial =
        reconstruct_cell(cells, right_cell, config);
    const Conservative1D left_average = sample_cell(cells, left_cell, config);
    const Conservative1D right_average = sample_cell(cells, right_cell, config);

    if (interface_is_troubled(left_average, right_average, config.gamma) ||
        !polynomial_is_admissible(left_polynomial, config.gamma) ||
        !polynomial_is_admissible(right_polynomial, config.gamma)) {
        return hll_flux(left_average, right_average, config.gamma);
    }

    Conservative1D flux_average{};
    try {
        for (const double tau : time_nodes) {
            Conservative1D left =
                predict_local_state(left_polynomial, 0.5, tau, dx, config.gamma);
            Conservative1D right =
                predict_local_state(right_polynomial, -0.5, tau, dx, config.gamma);

            if (!is_admissible(left, config.gamma) ||
                !is_admissible(right, config.gamma)) {
                return hll_flux(left_average, right_average, config.gamma);
            }
            flux_average += 0.5 * hll_flux(left, right, config.gamma);
        }
    } catch (const std::runtime_error&) {
        return hll_flux(left_average, right_average, config.gamma);
    }
    return flux_average;
}

}  // namespace

SimulationResult1D run_ader3_euler_1d(
    const RuntimeConfig& config,
    const std::vector<Conservative1D>& initial_state) {
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const double dx =
        (config.x_max - config.x_min) / static_cast<double>(config.nx);
    std::vector<Conservative1D> cells = initial_state;
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        double dt = compute_time_step(cells, config, dx);
        dt = std::min(dt, config.final_time - time);

        std::vector<Conservative1D> fluxes(config.nx + 1);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int iface = 0; iface <= config.nx; ++iface) {
            fluxes[iface] =
                time_averaged_interface_flux(cells, iface, config, dt, dx);
        }

        std::vector<Conservative1D> next(config.nx);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int i = 0; i < config.nx; ++i) {
            next[i] = cells[i] - (dt / dx) * (fluxes[i + 1] - fluxes[i]);
            if (!is_admissible(next[i], config.gamma)) {
                throw std::runtime_error(
                    "ADER-3 Euler update produced a non-admissible cell state.");
            }
        }
        cells = std::move(next);
        time += dt;
        ++step;
    }

    SimulationResult1D result;
    result.cells = std::move(cells);
    result.final_time = time;
    result.steps = step;
    result.x.resize(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] =
            config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }
    return result;
}

}  // namespace fp
