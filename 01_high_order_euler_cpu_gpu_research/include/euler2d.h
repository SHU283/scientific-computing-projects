#pragma once

#include <vector>

namespace fp {

// Conservative variables for the 2D Euler equations: U = (rho, rho*u, rho*v, E).
struct Conservative2D {
    double rho = 0.0;
    double mom_x = 0.0;
    double mom_y = 0.0;
    double energy = 0.0;
};

// Primitive variables: W = (rho, u, v, p).
struct Primitive2D {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double p = 0.0;
};

// Uniform Cartesian grid metadata.
struct Grid2D {
    int nx = 0;
    int ny = 0;
    double x_min = 0.0;
    double x_max = 1.0;
    double y_min = 0.0;
    double y_max = 1.0;
};

// Row-major cell-centered 2D field storage.
struct Field2D {
    int nx = 0;
    int ny = 0;
    std::vector<Conservative2D> cells;

    Conservative2D& operator()(int i, int j);
    const Conservative2D& operator()(int i, int j) const;
};

// Full output of a 2D solver run.
struct SimulationResult2D {
    Grid2D grid;
    Field2D field;
    double final_time = 0.0;
    int steps = 0;
};

// Conversions, physical fluxes, and Rusanov fluxes shared by 2D methods.
Primitive2D conservative_to_primitive(const Conservative2D& u, double gamma);
Conservative2D primitive_to_conservative(const Primitive2D& w, double gamma);
Conservative2D physical_flux_x(const Conservative2D& u, double gamma);
Conservative2D physical_flux_y(const Conservative2D& u, double gamma);
double sound_speed(const Primitive2D& w, double gamma);
double max_signal_speed_x(const Conservative2D& u, double gamma);
double max_signal_speed_y(const Conservative2D& u, double gamma);
Conservative2D operator+(const Conservative2D& a, const Conservative2D& b);
Conservative2D operator-(const Conservative2D& a, const Conservative2D& b);
Conservative2D operator*(double scalar, const Conservative2D& value);
Conservative2D operator*(const Conservative2D& value, double scalar);
Conservative2D operator/(const Conservative2D& value, double scalar);

// Limiter and Rusanov face fluxes used by the 2D finite-volume methods.
double minmod_2d(double a, double b);
Conservative2D rusanov_flux_x(const Conservative2D& left, const Conservative2D& right, double gamma);
Conservative2D rusanov_flux_y(const Conservative2D& bottom, const Conservative2D& top, double gamma);

// Uniform grid geometry helpers.
double dx(const Grid2D& grid);
double dy(const Grid2D& grid);
double x_center(const Grid2D& grid, int i);
double y_center(const Grid2D& grid, int j);

}  // namespace fp
