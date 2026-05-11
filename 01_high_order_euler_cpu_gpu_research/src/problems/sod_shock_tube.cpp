#include "problems.h"

#include <cmath>
#include <stdexcept>

namespace fp {

namespace {

double wrap_periodic(double x, double x_min, double x_max) {
    // Map a coordinate back into the periodic domain before evaluating exact
    // transported profiles.
    const double length = x_max - x_min;
    while (x < x_min) {
        x += length;
    }
    while (x >= x_max) {
        x -= length;
    }
    return x;
}

double smooth_density_average(double x_left,
                              double x_right,
                              double x_min,
                              double x_max,
                              double phase_shift) {
    // Return the cell average of the sinusoidal density profile. Using averages
    // instead of point samples makes convergence comparisons less grid-biased.
    const double length = x_max - x_min;
    constexpr double pi = 3.14159265358979323846;
    const double k = 2.0 * pi / length;

    const double integral =
        (-std::cos(k * (x_right - x_min - phase_shift)) + std::cos(k * (x_left - x_min - phase_shift))) / k;
    return 1.0 + 0.2 * integral / (x_right - x_left);
}

double gaussian_density_average(double x_left,
                                double x_right,
                                double center,
                                double sigma) {
    // Analytic integral of a Gaussian over one cell, expressed with erf.
    constexpr double sqrt_half = 0.70710678118654752440;
    constexpr double sqrt_pi_over_2 = 1.25331413731550025121;

    const double left_arg = (x_left - center) * sqrt_half / sigma;
    const double right_arg = (x_right - center) * sqrt_half / sigma;
    const double integral =
        sigma * sqrt_pi_over_2 * (std::erf(right_arg) - std::erf(left_arg));
    return 1.0 + 0.2 * integral / (x_right - x_left);
}

}  // namespace

std::vector<Conservative1D> make_initial_condition(const RuntimeConfig& config) {
    // Dispatcher used by the CLI and scripts so all 1D benchmarks share the same
    // initial-condition interface.
    if (config.problem == "sod") {
        return make_sod_initial_condition(config);
    }

    if (config.problem == "smooth_advection") {
        return make_smooth_advection_initial_condition(config);
    }

    if (config.problem == "gaussian_advection") {
        return make_gaussian_advection_initial_condition(config);
    }

    throw std::runtime_error("Unsupported problem: " + config.problem);
}

std::vector<Conservative1D> make_sod_initial_condition(const RuntimeConfig& config) {
    // Classic Sod shock tube: a left/right pressure and density jump, initially at rest.
    std::vector<Conservative1D> state(config.nx);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double x_mid = 0.5 * (config.x_min + config.x_max);

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;

        Primitive1D w;
        if (x < x_mid) {
            w = {1.0, 0.0, 1.0};
        } else {
            w = {0.125, 0.0, 0.1};
        }

        state[i] = primitive_to_conservative(w, config.gamma);
    }

    return state;
}

std::vector<Conservative1D> make_smooth_advection_initial_condition(const RuntimeConfig& config) {
    // Smooth periodic sine wave used mainly for convergence-order checks.
    std::vector<Conservative1D> state(config.nx);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double length = config.x_max - config.x_min;
    constexpr double pi = 3.14159265358979323846;
    const double phase_shift = 0.173 * length;

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;
        const double rho = 1.0 + 0.2 * std::sin(2.0 * pi * (x - config.x_min - phase_shift) / length);
        const Primitive1D w{rho, config.advection_velocity, 1.0};
        state[i] = primitive_to_conservative(w, config.gamma);
    }

    return state;
}

std::vector<Conservative1D> make_gaussian_advection_initial_condition(const RuntimeConfig& config) {
    // Smooth localized bump for testing numerical diffusion and periodic transport.
    std::vector<Conservative1D> state(config.nx);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double center = config.x_min + 0.317 * (config.x_max - config.x_min);
    const double sigma = 0.08 * (config.x_max - config.x_min);

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;
        const double r = (x - center) / sigma;
        const double rho = 1.0 + 0.2 * std::exp(-0.5 * r * r);
        const Primitive1D w{rho, config.advection_velocity, 1.0};
        state[i] = primitive_to_conservative(w, config.gamma);
    }

    return state;
}

SimulationResult1D make_exact_solution(const RuntimeConfig& config, double time) {
    // Exact density after pure advection. This is used only for smooth transport
    // benchmarks; discontinuous shock problems are evaluated by qualitative fields.
    SimulationResult1D exact;
    exact.x.resize(config.nx);
    exact.cells.resize(config.nx);
    exact.final_time = time;

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double length = config.x_max - config.x_min;
    constexpr double pi = 3.14159265358979323846;

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;
        exact.x[i] = x;
        const double x_left = x - 0.5 * dx;
        const double x_right = x + 0.5 * dx;

        if (config.problem == "smooth_advection") {
            const double phase_shift = 0.173 * length;
            
            const double shifted_left =
                wrap_periodic(x_left - config.advection_velocity * time, config.x_min, config.x_max);
            const double shifted_right =
                wrap_periodic(x_right - config.advection_velocity * time, config.x_min, config.x_max);

            double rho = 0.0;
            if (shifted_left <= shifted_right) {
                rho = smooth_density_average(shifted_left, shifted_right, config.x_min, config.x_max, phase_shift);
            } else {
                const double len_left = config.x_max - shifted_left;
                const double len_right = shifted_right - config.x_min;
                const double avg_left =
                    smooth_density_average(shifted_left, config.x_max, config.x_min, config.x_max, phase_shift);
                const double avg_right =
                    smooth_density_average(config.x_min, shifted_right, config.x_min, config.x_max, phase_shift);
                rho = (len_left * avg_left + len_right * avg_right) / (len_left + len_right);
            }

            const Primitive1D w{rho, config.advection_velocity, 1.0};
            exact.cells[i] = primitive_to_conservative(w, config.gamma);
        } else if (config.problem == "gaussian_advection") {
            const double shifted_center = wrap_periodic(
                config.x_min + 0.317 * (config.x_max - config.x_min) + config.advection_velocity * time,
                config.x_min,
                config.x_max);
            const double sigma = 0.08 * (config.x_max - config.x_min);

            double rho = 0.0;
            if (shifted_center - 4.0 * sigma >= config.x_min && shifted_center + 4.0 * sigma <= config.x_max) {
                rho = gaussian_density_average(x_left, x_right, shifted_center, sigma);
            } else {
                double center = shifted_center;
                if (x < config.x_min + 0.5 * length && shifted_center > config.x_min + 0.5 * length) {
                    center -= length;
                } else if (x > config.x_min + 0.5 * length && shifted_center < config.x_min + 0.5 * length) {
                    center += length;
                }
                rho = gaussian_density_average(x_left, x_right, center, sigma);
            }

            const Primitive1D w{rho, config.advection_velocity, 1.0};
            exact.cells[i] = primitive_to_conservative(w, config.gamma);
        } else {
            exact.cells[i] = Conservative1D{};
        }
    }

    return exact;
}

}  // namespace fp
