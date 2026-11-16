#include "problems2d.h"

#include <cmath>

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

double periodic_gaussian_average_1d(double left,
                                    double right,
                                    double center,
                                    double sigma,
                                    double domain_length) {
    const double scale = sigma * std::sqrt(0.5 * std::acos(-1.0));
    const double denominator = std::sqrt(2.0) * sigma;
    double integral = 0.0;

    for (int image = -2; image <= 2; ++image) {
        const double image_center = center + static_cast<double>(image) * domain_length;
        integral += scale *
                    (std::erf((right - image_center) / denominator) -
                     std::erf((left - image_center) / denominator));
    }
    return integral / (right - left);
}

double gaussian_cell_average(const Grid2D& grid,
                             int i,
                             int j,
                             double x0,
                             double y0,
                             double sigma) {
    const double hx = dx(grid);
    const double hy = dy(grid);
    const double x_left = grid.x_min + static_cast<double>(i) * hx;
    const double y_left = grid.y_min + static_cast<double>(j) * hy;
    const double average_x = periodic_gaussian_average_1d(
        x_left, x_left + hx, x0, sigma, grid.x_max - grid.x_min);
    const double average_y = periodic_gaussian_average_1d(
        y_left, y_left + hy, y0, sigma, grid.y_max - grid.y_min);
    return 1.0 + 0.2 * average_x * average_y;
}

}  // namespace

Field2D make_gaussian_advection_2d_initial_condition(const RuntimeConfig& config,
                                                     const Grid2D& grid) {
    Field2D field;
    field.nx = grid.nx;
    field.ny = grid.ny;
    field.cells.resize(static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny));

    const double x0 = grid.x_min + 0.31 * (grid.x_max - grid.x_min);
    const double y0 = grid.y_min + 0.43 * (grid.y_max - grid.y_min);
    const double sigma = 0.08 * (grid.x_max - grid.x_min);

    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const double rho = gaussian_cell_average(grid, i, j, x0, y0, sigma);
            const Primitive2D w{rho, config.advection_velocity, 0.37, 1.0};
            field(i, j) = primitive_to_conservative(w, config.gamma);
        }
    }

    return field;
}

Field2D make_gaussian_advection_2d_exact_solution(const RuntimeConfig& config,
                                                  const Grid2D& grid,
                                                  double time) {
    Field2D field;
    field.nx = grid.nx;
    field.ny = grid.ny;
    field.cells.resize(static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny));

    const double length_x = grid.x_max - grid.x_min;
    const double length_y = grid.y_max - grid.y_min;
    const double x0 = wrap_periodic(grid.x_min + 0.31 * length_x + config.advection_velocity * time,
                                    grid.x_min,
                                    grid.x_max);
    const double y0 = wrap_periodic(grid.y_min + 0.43 * length_y + 0.37 * time,
                                    grid.y_min,
                                    grid.y_max);
    const double sigma = 0.08 * length_x;

    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const double rho = gaussian_cell_average(grid, i, j, x0, y0, sigma);
            const Primitive2D w{rho, config.advection_velocity, 0.37, 1.0};
            field(i, j) = primitive_to_conservative(w, config.gamma);
        }
    }

    return field;
}

Field2D make_quadrant_riemann_2d_initial_condition(const RuntimeConfig& config,
                                                   const Grid2D& grid) {
    Field2D field;
    field.nx = grid.nx;
    field.ny = grid.ny;
    field.cells.resize(static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny));

    const double x_mid = 0.5 * (grid.x_min + grid.x_max);
    const double y_mid = 0.5 * (grid.y_min + grid.y_max);

    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const double x = x_center(grid, i);
            const double y = y_center(grid, j);

            Primitive2D w;
            if (x >= x_mid && y >= y_mid) {
                w = {1.5, 0.0, 0.0, 1.5};
            } else if (x < x_mid && y >= y_mid) {
                w = {0.5323, 1.206, 0.0, 0.3};
            } else if (x < x_mid && y < y_mid) {
                w = {0.138, 1.206, 1.206, 0.029};
            } else {
                w = {0.5323, 0.0, 1.206, 0.3};
            }

            field(i, j) = primitive_to_conservative(w, config.gamma);
        }
    }

    return field;
}

Field2D make_shock_bubble_2d_initial_condition(const RuntimeConfig& config,
                                               const Grid2D& grid) {
    Field2D field;
    field.nx = grid.nx;
    field.ny = grid.ny;
    field.cells.resize(static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny));

    const double bubble_x = grid.x_min + 0.35 * (grid.x_max - grid.x_min);
    const double bubble_y = grid.y_min + 0.5 * (grid.y_max - grid.y_min);
    const double bubble_radius = 0.12 * (grid.x_max - grid.x_min);

    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            const double x = x_center(grid, i);
            const double y = y_center(grid, j);

            Primitive2D w;
            if (x < grid.x_min + 0.15 * (grid.x_max - grid.x_min)) {
                w = {1.3764, 0.394, 0.0, 1.5698};
            } else {
                w = {1.0, 0.0, 0.0, 1.0};
            }

            const double dx_b = x - bubble_x;
            const double dy_b = y - bubble_y;
            const double r = std::sqrt(dx_b * dx_b + dy_b * dy_b);
            if (r <= bubble_radius) {
                w.rho = 0.1;
                w.p = 1.0;
                w.u = 0.0;
                w.v = 0.0;
            }

            field(i, j) = primitive_to_conservative(w, config.gamma);
        }
    }

    return field;
}

}  // namespace fp
