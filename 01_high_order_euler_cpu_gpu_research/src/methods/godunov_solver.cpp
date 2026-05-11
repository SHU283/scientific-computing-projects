#include "godunov_solver.h"

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

std::vector<Conservative1D> add_ghost_cells(const std::vector<Conservative1D>& interior) {
    // Legacy helper for outflow-style boundaries: copy edge states into one
    // ghost cell on each side.
    std::vector<Conservative1D> extended(interior.size() + 2);
    extended.front() = interior.front();
    extended.back() = interior.back();

    for (std::size_t i = 0; i < interior.size(); ++i) {
        extended[i + 1] = interior[i];
    }

    return extended;
}

std::vector<Conservative1D> add_ghost_cells(const std::vector<Conservative1D>& interior,
                                            const RuntimeConfig& config) {
    // One ghost cell on each side is sufficient for a first-order Godunov flux.
    std::vector<Conservative1D> extended(interior.size() + 2);

    for (std::size_t i = 0; i < interior.size(); ++i) {
        extended[i + 1] = interior[i];
    }

    if (config.boundary == "periodic") {
        extended.front() = interior.back();
        extended.back() = interior.front();
    } else {
        extended.front() = interior.front();
        extended.back() = interior.back();
    }

    return extended;
}

double compute_time_step(const std::vector<Conservative1D>& cells, const RuntimeConfig& config, double dx) {
    // CFL step size based on the largest local acoustic signal speed.
    double max_speed = 1.0e-12;
    for (std::size_t i = 0; i + 1 < cells.size(); ++i) {
        max_speed = std::max(max_speed, max_wave_speed(cells[i], cells[i + 1], config.gamma));
    }
    return config.cfl * dx / max_speed;
}

}  // namespace

SimulationResult1D run_godunov_1d(const RuntimeConfig& config,
                                  const std::vector<Conservative1D>& initial_state) {
    // First-order finite-volume baseline:
    // 1. extend with ghost cells,
    // 2. solve a Riemann problem at every cell face,
    // 3. update each cell by the conservative flux difference.
    if (initial_state.size() != static_cast<std::size_t>(config.nx)) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    std::vector<Conservative1D> cells = initial_state;

    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        const std::vector<Conservative1D> extended = add_ghost_cells(cells, config);
        double dt = compute_time_step(extended, config, dx);
        dt = std::min(dt, config.final_time - time);

        // Interface fluxes and cell updates are independent loops and can be OpenMP-parallelized.
        std::vector<Conservative1D> fluxes(config.nx + 1);
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int iface = 0; iface <= config.nx; ++iface) {
            fluxes[iface] = hll_flux(extended[iface], extended[iface + 1], config.gamma);
        }

        std::vector<Conservative1D> next = cells;
        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int i = 0; i < config.nx; ++i) {
            // U_i^{n+1} = U_i^n - dt/dx * (F_{i+1/2} - F_{i-1/2}).
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
