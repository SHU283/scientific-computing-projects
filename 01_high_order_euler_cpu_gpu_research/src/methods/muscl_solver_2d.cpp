#include "muscl_solver_2d.h"

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
    // Periodic benchmarks wrap indices; shock problems clamp to an outflow-style boundary.
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

Conservative2D limited_slope_x(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    // Directional minmod slopes provide second-order reconstruction while controlling shock oscillations.
    const Conservative2D& left = sample_cell(field, i - 1, j, config);
    const Conservative2D& center = sample_cell(field, i, j, config);
    const Conservative2D& right = sample_cell(field, i + 1, j, config);
    return {
        minmod_2d(center.rho - left.rho, right.rho - center.rho),
        minmod_2d(center.mom_x - left.mom_x, right.mom_x - center.mom_x),
        minmod_2d(center.mom_y - left.mom_y, right.mom_y - center.mom_y),
        minmod_2d(center.energy - left.energy, right.energy - center.energy),
    };
}

Conservative2D limited_slope_y(const Field2D& field, int i, int j, const RuntimeConfig& config) {
    const Conservative2D& bottom = sample_cell(field, i, j - 1, config);
    const Conservative2D& center = sample_cell(field, i, j, config);
    const Conservative2D& top = sample_cell(field, i, j + 1, config);
    return {
        minmod_2d(center.rho - bottom.rho, top.rho - center.rho),
        minmod_2d(center.mom_x - bottom.mom_x, top.mom_x - center.mom_x),
        minmod_2d(center.mom_y - bottom.mom_y, top.mom_y - center.mom_y),
        minmod_2d(center.energy - bottom.energy, top.energy - center.energy),
    };
}

Field2D add_scaled(const Field2D& base, const Field2D& delta, double scale) {
    // Utility for explicit Runge-Kutta stages: base + scale * RHS.
    Field2D result = base;
    for (std::size_t k = 0; k < base.cells.size(); ++k) {
        result.cells[k] = base.cells[k] + scale * delta.cells[k];
    }
    return result;
}

// Compute the semi-discrete finite-volume RHS:
// dU/dt = -dF/dx - dG/dy.
Field2D compute_rhs(const Field2D& field, const Grid2D& grid, const RuntimeConfig& config) {
    const double hx = dx(grid);
    const double hy = dy(grid);

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
            // Reconstruct left/right traces at each vertical face using x-slopes.
            const Conservative2D& left_cell = sample_cell(field, i - 1, j, config);
            const Conservative2D& right_cell = sample_cell(field, i, j, config);

            const Conservative2D slope_left = limited_slope_x(field, i - 1, j, config);
            const Conservative2D slope_right = limited_slope_x(field, i, j, config);

            const Conservative2D left_state = left_cell + 0.5 * slope_left;
            const Conservative2D right_state = right_cell - 0.5 * slope_right;
            fx_at(i, j) = rusanov_flux_x(left_state, right_state, config.gamma);
        }
    }

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j <= grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            // Reconstruct bottom/top traces at each horizontal face using y-slopes.
            const Conservative2D& bottom_cell = sample_cell(field, i, j - 1, config);
            const Conservative2D& top_cell = sample_cell(field, i, j, config);

            const Conservative2D slope_bottom = limited_slope_y(field, i, j - 1, config);
            const Conservative2D slope_top = limited_slope_y(field, i, j, config);

            const Conservative2D bottom_state = bottom_cell + 0.5 * slope_bottom;
            const Conservative2D top_state = top_cell - 0.5 * slope_top;
            fy_at(i, j) = rusanov_flux_y(bottom_state, top_state, config.gamma);
        }
    }

    Field2D rhs;
    rhs.nx = grid.nx;
    rhs.ny = grid.ny;
    rhs.cells.resize(field.cells.size());

    #pragma omp parallel for collapse(2) if(FP_HAS_OPENMP && config.backend == "openmp")
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            rhs(i, j) = (-1.0 / hx) * (fx_at(i + 1, j) - fx_at(i, j))
                      + (-1.0 / hy) * (fy_at(i, j + 1) - fy_at(i, j));
        }
    }

    return rhs;
}

}  // namespace

SimulationResult2D run_muscl_2d(const RuntimeConfig& config, const Grid2D& grid, const Field2D& initial_state) {
    // Second-order method-of-lines version of MUSCL in 2D. The spatial RHS uses
    // limited face states; SSP-RK2 provides the time update.
    Field2D field = initial_state;

    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        double dt = compute_time_step(field, grid, config);
        dt = std::min(dt, config.final_time - time);

        const Field2D rhs_1 = compute_rhs(field, grid, config);
        const Field2D stage_1 = add_scaled(field, rhs_1, dt);
        const Field2D rhs_2 = compute_rhs(stage_1, grid, config);

        #pragma omp parallel for if(FP_HAS_OPENMP && config.backend == "openmp")
        for (std::size_t k = 0; k < field.cells.size(); ++k) {
            // SSP-RK2 blend keeps the method stable for shock-capturing tests.
            field.cells[k] = 0.5 * field.cells[k] + 0.5 * (stage_1.cells[k] + dt * rhs_2.cells[k]);
        }

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
