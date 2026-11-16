#include "muscl_solver.h"

#include <algorithm>
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

// Add two transmissive ghost cells on each side for the MUSCL stencil.
std::vector<Conservative1D> add_ghost_cells(const std::vector<Conservative1D>& interior) {
    std::vector<Conservative1D> extended(interior.size() + 4);

    extended[0] = interior.front();
    extended[1] = interior.front();
    extended[extended.size() - 2] = interior.back();
    extended[extended.size() - 1] = interior.back();

    for (std::size_t i = 0; i < interior.size(); ++i) {
        extended[i + 2] = interior[i];
    }

    return extended;
}

// Add two ghost cells on each side using the configured boundary condition.
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

// CFL time step based on maximum interface wave speed.
double compute_time_step(const std::vector<Conservative1D>& cells, const RuntimeConfig& config, double dx) {
    double max_speed = 1.0e-12;
    for (std::size_t i = 0; i + 1 < cells.size(); ++i) {
        max_speed = std::max(max_speed, max_wave_speed(cells[i], cells[i + 1], config.gamma));
    }
    return config.cfl * dx / max_speed;
}

// Component-wise minmod slope used by MUSCL reconstruction.
Conservative1D limited_slope(const Conservative1D& left_diff, const Conservative1D& right_diff) {
    return {
        minmod(left_diff.rho, right_diff.rho),
        minmod(left_diff.mom, right_diff.mom),
        minmod(left_diff.energy, right_diff.energy),
    };
}

}  // namespace

SimulationResult1D run_muscl_hancock_1d(const RuntimeConfig& config,
                                        const std::vector<Conservative1D>& initial_state) {
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    std::vector<Conservative1D> cells = initial_state;

    // Initialize time and step counters.
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        // Extended state supplies the two-cell stencil near the boundaries.
        const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);
        double dt = compute_time_step(extended, config, dx);
        // Clip the last step to land exactly on final_time.
        dt = std::min(dt, config.final_time - time);

        // Compute limited slopes for every cell in the extended grid.
        std::vector<Conservative1D> slopes(config.nx + 4);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int i = 1; i < config.nx + 3; ++i) {
            const Conservative1D left_diff = extended[i] - extended[i - 1];
            const Conservative1D right_diff = extended[i + 1] - extended[i];
            slopes[i] = limited_slope(left_diff, right_diff);
        }

        std::vector<Conservative1D> left_states(config.nx + 1);
        std::vector<Conservative1D> right_states(config.nx + 1);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int iface = 0; iface <= config.nx; ++iface) {
            // Interface maps to adjacent extended-grid cell indices.
            const int cell_left = iface + 1;
            const int cell_right = iface + 2;

            const Conservative1D u_left = extended[cell_left];
            const Conservative1D u_right = extended[cell_right];

            // Physical fluxes at each side of each cell provide the Hancock
            // half-step predictor.
            const Conservative1D left_flux_plus = physical_flux(u_left + 0.5 * slopes[cell_left], config.gamma);
            const Conservative1D left_flux_minus = physical_flux(u_left - 0.5 * slopes[cell_left], config.gamma);
            const Conservative1D right_flux_plus = physical_flux(u_right + 0.5 * slopes[cell_right], config.gamma);
            const Conservative1D right_flux_minus = physical_flux(u_right - 0.5 * slopes[cell_right], config.gamma);

            left_states[iface] =
                (u_left + 0.5 * slopes[cell_left]) - 0.5 * (dt / dx) * (left_flux_plus - left_flux_minus);
            right_states[iface] =
                (u_right - 0.5 * slopes[cell_right]) - 0.5 * (dt / dx) * (right_flux_plus - right_flux_minus);
        }

        // HLL flux from the predicted left/right interface states.
        std::vector<Conservative1D> fluxes(config.nx + 1);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int iface = 0; iface <= config.nx; ++iface) {
            fluxes[iface] = hll_flux(left_states[iface], right_states[iface], config.gamma);
        }

        std::vector<Conservative1D> next = cells;
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int i = 0; i < config.nx; ++i) {
            next[i] = cells[i] - (dt / dx) * (fluxes[i + 1] - fluxes[i]);
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

    #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }

    return result;
}

}  // namespace fp

