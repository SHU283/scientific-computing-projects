#include "godunov_solver_2d.h"

#include <algorithm>
#include <vector>

#ifdef FP_HAS_OPENMP
#include <omp.h>
#endif

#ifndef FP_HAS_OPENMP
#define FP_HAS_OPENMP 0
#endif

namespace fp {

namespace {

int wrap_index(int idx, int size) {
    while (idx < 0) {
        idx += size;
    }
    while (idx >= size) {
        idx -= size;
    }
    return idx;
}
 
const Conservative2D& sample_cell(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    if (config.boundary == "periodic") {
        return field(wrap_index(i, field.nx), wrap_index(j, field.ny));
    }

    i = std::max(0, std::min(i, field.nx - 1));
    j = std::max(0, std::min(j, field.ny - 1));
    return field(i, j);
}

double compute_time_step(const Field2D& field, const Grid2D& grid, const RuntimeConfig& config) {
    double max_speed = 1.0e-12;
    for (int j = 0; j < field.ny; ++j) {
        for (int i = 0; i < field.nx; ++i) {
            max_speed = std::max(max_speed, max_signal_speed_x(field(i, j), config.gamma));
            max_speed = std::max(max_speed, max_signal_speed_y(field(i, j), config.gamma));
        }
    }
    return config.cfl * std::min(dx(grid), dy(grid)) / max_speed;
}

}  // namespace

SimulationResult2D run_godunov_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state) {
    Field2D field = initial_state;

    double time = 0.0;
    int step = 0;
    const double hx = dx(grid);
    const double hy = dy(grid);

    while (time < config.final_time && step < config.max_steps) {
        double dt = compute_time_step(field, grid, config);
        dt = std::min(dt, config.final_time - time);

        std::vector<Conservative2D> fx(static_cast<std::size_t>((grid.nx + 1) * grid.ny));
        std::vector<Conservative2D> fy(static_cast<std::size_t>(grid.nx * (grid.ny + 1)));

        auto fx_at = [&](int i, int j) -> Conservative2D& {
            return fx[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx + 1) + static_cast<std::size_t>(i)];
        };
        auto fy_at = [&](int i, int j) -> Conservative2D& {
            return fy[static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx) + static_cast<std::size_t>(i)];
        };

        #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i <= grid.nx; ++i) {
                const Conservative2D& left = sample_cell(field, i - 1, j, config);
                const Conservative2D& right = sample_cell(field, i, j, config);
                fx_at(i, j) = rusanov_flux_x(left, right, config.gamma);
            }
        }

        #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int j = 0; j <= grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                const Conservative2D& bottom = sample_cell(field, i, j - 1, config);
                const Conservative2D& top = sample_cell(field, i, j, config);
                fy_at(i, j) = rusanov_flux_y(bottom, top, config.gamma);
            }
        }

        Field2D next = field;
        #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                const Conservative2D update_x = (dt / hx) * (fx_at(i + 1, j) - fx_at(i, j));
                const Conservative2D update_y = (dt / hy) * (fy_at(i, j + 1) - fy_at(i, j));
                next(i, j) = field(i, j) - update_x - update_y;
            }
        }

        field = std::move(next);
        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field = std::move(field);
    result.final_time = time;
    result.steps = step;
    return result;
}

}  // namespace fp
