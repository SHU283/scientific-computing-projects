#include "ader3_euler_solver_2d.h"

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
constexpr double kTroubledJumpThreshold = 0.15;
constexpr int kMaximumStepRetries = 8;

int sample_index(int index, int size, const RuntimeConfig& config) {
    if (config.boundary == "periodic") {
        index %= size;
        return index < 0 ? index + size : index;
    }
    return std::clamp(index, 0, size - 1);
}

Conservative2D sample_cell(const Field2D& field,
                           int i,
                           int j,
                           const RuntimeConfig& config) {
    return field(
        sample_index(i, field.nx, config),
        sample_index(j, field.ny, config));
}

bool is_admissible(const Conservative2D& state, double gamma) {
    if (!(state.rho > 0.0) || !std::isfinite(state.rho) ||
        !std::isfinite(state.mom_x) || !std::isfinite(state.mom_y) ||
        !std::isfinite(state.energy)) {
        return false;
    }
    const double kinetic =
        0.5 * (state.mom_x * state.mom_x + state.mom_y * state.mom_y) /
        state.rho;
    return (gamma - 1.0) * (state.energy - kinetic) > 0.0;
}

struct LocalPolynomial2D {
    Conservative2D constant;
    Conservative2D slope_x;
    Conservative2D slope_y;
    Conservative2D curvature_x;
    Conservative2D curvature_xy;
    Conservative2D curvature_y;

    Conservative2D evaluate(double xi, double eta) const {
        return constant +
               xi * slope_x +
               eta * slope_y +
               (0.5 * xi * xi) * curvature_x +
               (xi * eta) * curvature_xy +
               (0.5 * eta * eta) * curvature_y;
    }
};

LocalPolynomial2D reconstruct_cell(const Field2D& field,
                                   int i,
                                   int j,
                                   const RuntimeConfig& config) {
    const Conservative2D center = sample_cell(field, i, j, config);
    const Conservative2D left = sample_cell(field, i - 1, j, config);
    const Conservative2D right = sample_cell(field, i + 1, j, config);
    const Conservative2D bottom = sample_cell(field, i, j - 1, config);
    const Conservative2D top = sample_cell(field, i, j + 1, config);

    const Conservative2D curvature_x = right - 2.0 * center + left;
    const Conservative2D curvature_y = top - 2.0 * center + bottom;
    const Conservative2D curvature_xy =
        0.25 * (
            sample_cell(field, i + 1, j + 1, config) -
            sample_cell(field, i + 1, j - 1, config) -
            sample_cell(field, i - 1, j + 1, config) +
            sample_cell(field, i - 1, j - 1, config));

    return {
        center - (1.0 / 24.0) * (curvature_x + curvature_y),
        0.5 * (right - left),
        0.5 * (top - bottom),
        curvature_x,
        curvature_xy,
        curvature_y,
    };
}

bool polynomial_is_admissible(const LocalPolynomial2D& polynomial,
                              double gamma) {
    constexpr std::array<double, 5> points = {
        -0.5 - 2.0 * kReferenceDerivativeStep,
        -0.5,
        0.0,
        0.5,
        0.5 + 2.0 * kReferenceDerivativeStep,
    };
    for (const double eta : points) {
        for (const double xi : points) {
            if (!is_admissible(polynomial.evaluate(xi, eta), gamma)) {
                return false;
            }
        }
    }
    return true;
}

double max_component_scale(const Conservative2D& state) {
    return std::max({
        1.0,
        std::abs(state.rho),
        std::abs(state.mom_x),
        std::abs(state.mom_y),
        std::abs(state.energy),
    });
}

template <typename FluxFunction>
Conservative2D flux_jacobian_vector_product(
    const Conservative2D& state,
    const Conservative2D& direction,
    double gamma,
    FluxFunction flux) {
    const double state_scale = max_component_scale(state);
    const double direction_scale = max_component_scale(direction);
    double epsilon =
        kDirectionalDifferenceScale * state_scale / direction_scale;

    Conservative2D plus = state + epsilon * direction;
    Conservative2D minus = state - epsilon * direction;
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
            "ADER-3 2D predictor could not form an admissible Jacobian perturbation.");
    }
    return (flux(plus, gamma) - flux(minus, gamma)) / (2.0 * epsilon);
}

Conservative2D spatial_time_derivative(
    const LocalPolynomial2D& polynomial,
    double xi,
    double eta,
    double cell_dx,
    double cell_dy,
    double gamma) {
    const double h = kReferenceDerivativeStep;
    const Conservative2D flux_x_plus =
        physical_flux_x(polynomial.evaluate(xi + h, eta), gamma);
    const Conservative2D flux_x_minus =
        physical_flux_x(polynomial.evaluate(xi - h, eta), gamma);
    const Conservative2D flux_y_plus =
        physical_flux_y(polynomial.evaluate(xi, eta + h), gamma);
    const Conservative2D flux_y_minus =
        physical_flux_y(polynomial.evaluate(xi, eta - h), gamma);

    return (-1.0 / (2.0 * h * cell_dx)) *
               (flux_x_plus - flux_x_minus) +
           (-1.0 / (2.0 * h * cell_dy)) *
               (flux_y_plus - flux_y_minus);
}

Conservative2D second_time_derivative(
    const LocalPolynomial2D& polynomial,
    double xi,
    double eta,
    double cell_dx,
    double cell_dy,
    double gamma) {
    const double h = kReferenceDerivativeStep;

    const Conservative2D state_x_plus = polynomial.evaluate(xi + h, eta);
    const Conservative2D state_x_minus = polynomial.evaluate(xi - h, eta);
    const Conservative2D ut_x_plus = spatial_time_derivative(
        polynomial, xi + h, eta, cell_dx, cell_dy, gamma);
    const Conservative2D ut_x_minus = spatial_time_derivative(
        polynomial, xi - h, eta, cell_dx, cell_dy, gamma);
    const Conservative2D a_ut_plus = flux_jacobian_vector_product(
        state_x_plus, ut_x_plus, gamma, physical_flux_x);
    const Conservative2D a_ut_minus = flux_jacobian_vector_product(
        state_x_minus, ut_x_minus, gamma, physical_flux_x);

    const Conservative2D state_y_plus = polynomial.evaluate(xi, eta + h);
    const Conservative2D state_y_minus = polynomial.evaluate(xi, eta - h);
    const Conservative2D ut_y_plus = spatial_time_derivative(
        polynomial, xi, eta + h, cell_dx, cell_dy, gamma);
    const Conservative2D ut_y_minus = spatial_time_derivative(
        polynomial, xi, eta - h, cell_dx, cell_dy, gamma);
    const Conservative2D b_ut_plus = flux_jacobian_vector_product(
        state_y_plus, ut_y_plus, gamma, physical_flux_y);
    const Conservative2D b_ut_minus = flux_jacobian_vector_product(
        state_y_minus, ut_y_minus, gamma, physical_flux_y);

    return (-1.0 / (2.0 * h * cell_dx)) *
               (a_ut_plus - a_ut_minus) +
           (-1.0 / (2.0 * h * cell_dy)) *
               (b_ut_plus - b_ut_minus);
}

Conservative2D predict_local_state(
    const LocalPolynomial2D& polynomial,
    double xi,
    double eta,
    double tau,
    double cell_dx,
    double cell_dy,
    double gamma) {
    const Conservative2D state = polynomial.evaluate(xi, eta);
    const Conservative2D ut = spatial_time_derivative(
        polynomial, xi, eta, cell_dx, cell_dy, gamma);
    const Conservative2D utt = second_time_derivative(
        polynomial, xi, eta, cell_dx, cell_dy, gamma);
    return state + tau * ut + (0.5 * tau * tau) * utt;
}

bool interface_is_troubled(const Conservative2D& first,
                           const Conservative2D& second,
                           double gamma,
                           bool normal_is_x) {
    const Primitive2D left = conservative_to_primitive(first, gamma);
    const Primitive2D right = conservative_to_primitive(second, gamma);
    const double rho_jump =
        std::abs(right.rho - left.rho) /
        std::max({left.rho, right.rho, 1.0e-12});
    const double pressure_jump =
        std::abs(right.p - left.p) /
        std::max({left.p, right.p, 1.0e-12});
    const double sound_scale =
        std::max(sound_speed(left, gamma), sound_speed(right, gamma));
    const double normal_velocity_left = normal_is_x ? left.u : left.v;
    const double normal_velocity_right = normal_is_x ? right.u : right.v;
    const double velocity_jump =
        std::abs(normal_velocity_right - normal_velocity_left) /
        std::max(sound_scale, 1.0e-12);
    return std::max({rho_jump, pressure_jump, velocity_jump}) >
           kTroubledJumpThreshold;
}

double compute_time_step(const Field2D& field,
                         const Grid2D& grid,
                         const RuntimeConfig& config) {
    double max_x = 1.0e-12;
    double max_y = 1.0e-12;
    #pragma omp parallel for collapse(2) reduction(max:max_x,max_y) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < field.ny; ++j) {
        for (int i = 0; i < field.nx; ++i) {
            max_x = std::max(
                max_x, max_signal_speed_x(field(i, j), config.gamma));
            max_y = std::max(
                max_y, max_signal_speed_y(field(i, j), config.gamma));
        }
    }
    return config.cfl / (max_x / dx(grid) + max_y / dy(grid));
}

Conservative2D time_averaged_x_flux(
    const Field2D& field,
    int iface,
    int j,
    const RuntimeConfig& config,
    double dt,
    double cell_dx,
    double cell_dy) {
    constexpr double gauss = 0.28867513459481288225;
    constexpr std::array<double, 2> face_nodes = {-gauss, gauss};
    const std::array<double, 2> time_nodes = {
        dt * (0.5 - gauss),
        dt * (0.5 + gauss),
    };

    const int left_i = iface - 1;
    const int right_i = iface;
    const Conservative2D left_average =
        sample_cell(field, left_i, j, config);
    const Conservative2D right_average =
        sample_cell(field, right_i, j, config);
    const LocalPolynomial2D left_polynomial =
        reconstruct_cell(field, left_i, j, config);
    const LocalPolynomial2D right_polynomial =
        reconstruct_cell(field, right_i, j, config);

    if (interface_is_troubled(
            left_average, right_average, config.gamma, true) ||
        !polynomial_is_admissible(left_polynomial, config.gamma) ||
        !polynomial_is_admissible(right_polynomial, config.gamma)) {
        return rusanov_flux_x(
            left_average, right_average, config.gamma);
    }

    Conservative2D flux_average{};
    try {
        for (const double eta : face_nodes) {
            for (const double tau : time_nodes) {
                const Conservative2D left = predict_local_state(
                    left_polynomial, 0.5, eta, tau,
                    cell_dx, cell_dy, config.gamma);
                const Conservative2D right = predict_local_state(
                    right_polynomial, -0.5, eta, tau,
                    cell_dx, cell_dy, config.gamma);
                if (!is_admissible(left, config.gamma) ||
                    !is_admissible(right, config.gamma)) {
                    return rusanov_flux_x(
                        left_average, right_average, config.gamma);
                }
                flux_average = flux_average +
                    0.25 * rusanov_flux_x(left, right, config.gamma);
            }
        }
    } catch (const std::runtime_error&) {
        return rusanov_flux_x(left_average, right_average, config.gamma);
    }
    return flux_average;
}

Conservative2D time_averaged_y_flux(
    const Field2D& field,
    int i,
    int jface,
    const RuntimeConfig& config,
    double dt,
    double cell_dx,
    double cell_dy) {
    constexpr double gauss = 0.28867513459481288225;
    constexpr std::array<double, 2> face_nodes = {-gauss, gauss};
    const std::array<double, 2> time_nodes = {
        dt * (0.5 - gauss),
        dt * (0.5 + gauss),
    };

    const int bottom_j = jface - 1;
    const int top_j = jface;
    const Conservative2D bottom_average =
        sample_cell(field, i, bottom_j, config);
    const Conservative2D top_average =
        sample_cell(field, i, top_j, config);
    const LocalPolynomial2D bottom_polynomial =
        reconstruct_cell(field, i, bottom_j, config);
    const LocalPolynomial2D top_polynomial =
        reconstruct_cell(field, i, top_j, config);

    if (interface_is_troubled(
            bottom_average, top_average, config.gamma, false) ||
        !polynomial_is_admissible(bottom_polynomial, config.gamma) ||
        !polynomial_is_admissible(top_polynomial, config.gamma)) {
        return rusanov_flux_y(
            bottom_average, top_average, config.gamma);
    }

    Conservative2D flux_average{};
    try {
        for (const double xi : face_nodes) {
            for (const double tau : time_nodes) {
                const Conservative2D bottom = predict_local_state(
                    bottom_polynomial, xi, 0.5, tau,
                    cell_dx, cell_dy, config.gamma);
                const Conservative2D top = predict_local_state(
                    top_polynomial, xi, -0.5, tau,
                    cell_dx, cell_dy, config.gamma);
                if (!is_admissible(bottom, config.gamma) ||
                    !is_admissible(top, config.gamma)) {
                    return rusanov_flux_y(
                        bottom_average, top_average, config.gamma);
                }
                flux_average = flux_average +
                    0.25 * rusanov_flux_y(bottom, top, config.gamma);
            }
        }
    } catch (const std::runtime_error&) {
        return rusanov_flux_y(
            bottom_average, top_average, config.gamma);
    }
    return flux_average;
}

bool attempt_step(const Field2D& field,
                  Field2D& next,
                  const Grid2D& grid,
                  const RuntimeConfig& config,
                  double dt) {
    std::vector<Conservative2D> flux_x(
        static_cast<std::size_t>(grid.nx + 1) *
        static_cast<std::size_t>(grid.ny));
    std::vector<Conservative2D> flux_y(
        static_cast<std::size_t>(grid.nx) *
        static_cast<std::size_t>(grid.ny + 1));

    auto x_flux = [&](int iface, int j) -> Conservative2D& {
        return flux_x[
            static_cast<std::size_t>(j) *
                static_cast<std::size_t>(grid.nx + 1) +
            static_cast<std::size_t>(iface)];
    };
    auto y_flux = [&](int i, int jface) -> Conservative2D& {
        return flux_y[
            static_cast<std::size_t>(jface) *
                static_cast<std::size_t>(grid.nx) +
            static_cast<std::size_t>(i)];
    };

    const double cell_dx = dx(grid);
    const double cell_dy = dy(grid);
    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int iface = 0; iface <= grid.nx; ++iface) {
            x_flux(iface, j) = time_averaged_x_flux(
                field, iface, j, config, dt, cell_dx, cell_dy);
        }
    }
    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int jface = 0; jface <= grid.ny; ++jface) {
        for (int i = 0; i < grid.nx; ++i) {
            y_flux(i, jface) = time_averaged_y_flux(
                field, i, jface, config, dt, cell_dx, cell_dy);
        }
    }

    next = field;
    int invalid_cells = 0;
    #pragma omp parallel for collapse(2) reduction(+:invalid_cells) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            next(i, j) =
                field(i, j) -
                (dt / cell_dx) *
                    (x_flux(i + 1, j) - x_flux(i, j)) -
                (dt / cell_dy) *
                    (y_flux(i, j + 1) - y_flux(i, j));
            if (!is_admissible(next(i, j), config.gamma)) {
                ++invalid_cells;
            }
        }
    }
    return invalid_cells == 0;
}

}  // namespace

SimulationResult2D run_ader3_euler_2d(
    const RuntimeConfig& config,
    const Grid2D& grid,
    const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error(
            "Initial 2D field does not match the configured grid.");
    }

    Field2D field = initial_state;
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        double dt = std::min(
            compute_time_step(field, grid, config),
            config.final_time - time);

        Field2D next;
        bool accepted = false;
        for (int retry = 0; retry <= kMaximumStepRetries; ++retry) {
            if (attempt_step(field, next, grid, config, dt)) {
                accepted = true;
                break;
            }
            dt *= 0.5;
        }
        if (!accepted) {
            throw std::runtime_error(
                "ADER-3 2D could not produce an admissible update after reducing dt.");
        }

        field = std::move(next);
        time += dt;
        ++step;
    }

    if (step >= config.max_steps && time < config.final_time) {
        throw std::runtime_error(
            "ADER-3 2D reached max_steps before final_time.");
    }
    return {grid, std::move(field), time, step};
}

}  // namespace fp
