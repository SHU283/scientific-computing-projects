#include "problems.h"

#include <cmath>
#include <stdexcept>

namespace fp {

namespace {

double wrap_periodic(double x, double x_min, double x_max) {
    const double length = x_max - x_min;
    while (x < x_min) {
        x += length;
    }
    while (x >= x_max) {
        x -= length;
    }
    return x;
}

// -----------------------------------------------------------------------------
// Smooth sinusoidal density: exact finite-volume cell average
// -----------------------------------------------------------------------------
//
// Continuous density profile:
//
//     rho(x) = 1 + 0.2 * sin(k * (x - x_min - phase_shift))
//
// where
//
//     L = x_max - x_min
//     k = 2*pi / L
//
// Here L is the length of the full computational domain, so the sine wave
// completes one full period over [x_min, x_max]. The parameter phase_shift
// represents a spatial translation of the wave, for example phase_shift = a*t
// for constant-speed advection.
//
// In a finite-volume method, the stored value in each cell is the cell average,
// not the point value:
//
//     rho_bar_i = (1 / Delta x_i) * integral_{x_left}^{x_right} rho(x) dx
//
// where
//
//     Delta x_i = x_right - x_left
//
// Substituting the sine density:
//
//     rho_bar_i
//       = 1
//         + 0.2 / (x_right - x_left)
//           * integral_{x_left}^{x_right}
//               sin(k * (x - x_min - phase_shift)) dx
//
// Since
//
//     integral sin(k * (x - x_min - phase_shift)) dx
//       = -cos(k * (x - x_min - phase_shift)) / k
//
// the exact integral over the cell is
//
//     integral_{x_left}^{x_right} sin(k * (x - x_min - phase_shift)) dx
//       = [ -cos(k * (x - x_min - phase_shift)) / k ]_{x_left}^{x_right}
//
//       = ( -cos(k * (x_right - x_min - phase_shift))
//           + cos(k * (x_left  - x_min - phase_shift)) ) / k
//
// Therefore the exact cell-average density is
//
//     rho_bar_i
//       = 1
//         + 0.2
//           * ( -cos(k * (x_right - x_min - phase_shift))
//               + cos(k * (x_left  - x_min - phase_shift)) )
//           / ( k * (x_right - x_left) )
//
// This is what smooth_density_average(...) computes.
// -----------------------------------------------------------------------------
double smooth_density_average(double x_left,
                              double x_right,
                              double x_min,
                              double x_max,
                              double phase_shift) {
    const double length = x_max - x_min;
    constexpr double pi = 3.14159265358979323846;
    const double k = 2.0 * pi / length;

    const double integral =
        (-std::cos(k * (x_right - x_min - phase_shift)) + std::cos(k * (x_left - x_min - phase_shift))) / k;
    return 1.0 + 0.2 * integral / (x_right - x_left);
}

// -----------------------------------------------------------------------------
// Gaussian density perturbation: exact finite-volume cell average
// -----------------------------------------------------------------------------
//
// Continuous density profile:
//
//     rho(x) = 1 + 0.2 * exp( - (x - center)^2 / (2*sigma^2) )
//
// where
//
//     center = centre position of the Gaussian bump
//     sigma  = width parameter / standard deviation of the Gaussian
//
// In a finite-volume method, the stored density in each cell is the cell average:
//
//     rho_bar_i = (1 / Delta x_i) * integral_{x_left}^{x_right} rho(x) dx
//
// where
//
//     Delta x_i = x_right - x_left
//
// Substituting the Gaussian density:
//
//     rho_bar_i
//       = 1
//         + 0.2 / (x_right - x_left)
//           * integral_{x_left}^{x_right}
//               exp( - (x - center)^2 / (2*sigma^2) ) dx
//
// The Gaussian integral is expressed using the error function erf:
//
//     erf(z) = (2 / sqrt(pi)) * integral_0^z exp(-t^2) dt
//
// To derive the cell integral, use the change of variables:
//
//     z = (x - center) / (sqrt(2) * sigma)
//
// Then
//
//     dx = sqrt(2) * sigma dz
//
// and
//
//     exp( - (x - center)^2 / (2*sigma^2) ) = exp(-z^2)
//
// Therefore
//
//     integral exp( - (x - center)^2 / (2*sigma^2) ) dx
//       = sigma * sqrt(pi / 2)
//         * erf( (x - center) / (sqrt(2) * sigma) )
//
// Hence the exact integral over [x_left, x_right] is
//
//     integral_{x_left}^{x_right}
//         exp( - (x - center)^2 / (2*sigma^2) ) dx
//
//       = sigma * sqrt(pi / 2)
//         * [
//             erf( (x_right - center) / (sqrt(2) * sigma) )
//             -
//             erf( (x_left  - center) / (sqrt(2) * sigma) )
//           ]
//
// Therefore the exact cell-average density is
//
//     rho_bar_i
//       = 1
//         + 0.2
//           * sigma * sqrt(pi / 2)
//           * [
//               erf( (x_right - center) / (sqrt(2) * sigma) )
//               -
//               erf( (x_left  - center) / (sqrt(2) * sigma) )
//             ]
//           / (x_right - x_left)
//
// In the code:
//
//     sqrt_half      = 1 / sqrt(2)
//     sqrt_pi_over_2 = sqrt(pi / 2)
//
// This is what gaussian_density_average(...) computes.
// -----------------------------------------------------------------------------
double gaussian_density_average(double x_left,
                                double x_right,
                                double center,
                                double sigma) {
    constexpr double sqrt_half = 0.70710678118654752440;
    constexpr double sqrt_pi_over_2 = 1.25331413731550025121;

    const double left_arg = (x_left - center) * sqrt_half / sigma;
    const double right_arg = (x_right - center) * sqrt_half / sigma;
    const double integral =
        sigma * sqrt_pi_over_2 * (std::erf(right_arg) - std::erf(left_arg));
    return 1.0 + 0.2 * integral / (x_right - x_left);
}

double periodic_gaussian_density_average(double x_left,
                                         double x_right,
                                         double center,
                                         double sigma,
                                         double period) {
    double rho = 1.0;
    for (const int image : {-1, 0, 1}) {
        rho += gaussian_density_average(
                   x_left,
                   x_right,
                   center + static_cast<double>(image) * period,
                   sigma) -
               1.0;
    }
    return rho;
}

}  // namespace

std::vector<Conservative1D> make_initial_condition(const RuntimeConfig& config) {
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
    std::vector<Conservative1D> state(config.nx);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double length = config.x_max - config.x_min;
    const double phase_shift = 0.173 * length;

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;
        const double x_left = x - 0.5 * dx;
        const double x_right = x + 0.5 * dx;
        const double rho =
            smooth_density_average(x_left, x_right, config.x_min, config.x_max, phase_shift);
        const Primitive1D w{rho, config.advection_velocity, 1.0};
        state[i] = primitive_to_conservative(w, config.gamma);
    }

    return state;
}

std::vector<Conservative1D> make_gaussian_advection_initial_condition(const RuntimeConfig& config) {
    std::vector<Conservative1D> state(config.nx);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const double center = config.x_min + 0.317 * (config.x_max - config.x_min);
    const double sigma = 0.08 * (config.x_max - config.x_min);

    for (int i = 0; i < config.nx; ++i) {
        const double x = config.x_min + (static_cast<double>(i) + 0.5) * dx;
        const double x_left = x - 0.5 * dx;
        const double x_right = x + 0.5 * dx;
        const double rho = periodic_gaussian_density_average(
            x_left, x_right, center, sigma, config.x_max - config.x_min);
        const Primitive1D w{rho, config.advection_velocity, 1.0};
        state[i] = primitive_to_conservative(w, config.gamma);
    }

    return state;
}

SimulationResult1D make_exact_solution(const RuntimeConfig& config, double time) {
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

            const double rho = periodic_gaussian_density_average(
                x_left, x_right, shifted_center, sigma, length);

            const Primitive1D w{rho, config.advection_velocity, 1.0};
            exact.cells[i] = primitive_to_conservative(w, config.gamma);
        } else {
            exact.cells[i] = Conservative1D{};
        }
    }

    return exact;
}

}  // namespace fp
