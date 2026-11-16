#include "euler2d.h"

#include <cmath>
#include <stdexcept>

namespace fp {

Conservative2D& Field2D::operator()(int i, int j) {
    return cells[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i)];
}

const Conservative2D& Field2D::operator()(int i, int j) const {
    return cells[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i)];
}

Primitive2D conservative_to_primitive(const Conservative2D& u, double gamma) {
    if (u.rho <= 0.0) {
        throw std::runtime_error("2D density must stay positive.");
    }

    const double vel_x = u.mom_x / u.rho;
    const double vel_y = u.mom_y / u.rho;
    const double kinetic = 0.5 * u.rho * (vel_x * vel_x + vel_y * vel_y);
    const double pressure = (gamma - 1.0) * (u.energy - kinetic);

    if (pressure <= 0.0) {
        throw std::runtime_error("2D pressure must stay positive.");
    }

    return {u.rho, vel_x, vel_y, pressure};
}

Conservative2D operator+(const Conservative2D& a, const Conservative2D& b) {
    return {a.rho + b.rho, a.mom_x + b.mom_x, a.mom_y + b.mom_y, a.energy + b.energy};
}

Conservative2D operator-(const Conservative2D& a, const Conservative2D& b) {
    return {a.rho - b.rho, a.mom_x - b.mom_x, a.mom_y - b.mom_y, a.energy - b.energy};
}

Conservative2D operator*(double scalar, const Conservative2D& value) {
    return {scalar * value.rho, scalar * value.mom_x, scalar * value.mom_y, scalar * value.energy};
}

Conservative2D operator*(const Conservative2D& value, double scalar) {
    return scalar * value;
}

Conservative2D operator/(const Conservative2D& value, double scalar) {
    return {value.rho / scalar, value.mom_x / scalar, value.mom_y / scalar, value.energy / scalar};
}

double minmod_2d(double a, double b) {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return (std::abs(a) < std::abs(b)) ? a : b;
}

Conservative2D primitive_to_conservative(const Primitive2D& w, double gamma) {
    const double mom_x = w.rho * w.u;
    const double mom_y = w.rho * w.v;
    const double internal = w.p / (gamma - 1.0);
    const double kinetic = 0.5 * w.rho * (w.u * w.u + w.v * w.v);
    return {w.rho, mom_x, mom_y, internal + kinetic};
}

Conservative2D physical_flux_x(const Conservative2D& u, double gamma) {
    const Primitive2D w = conservative_to_primitive(u, gamma);
    return {
        u.mom_x,
        u.mom_x * w.u + w.p,
        u.mom_y * w.u,
        (u.energy + w.p) * w.u,
    };
}

Conservative2D physical_flux_y(const Conservative2D& u, double gamma) {
    const Primitive2D w = conservative_to_primitive(u, gamma);
    return {
        u.mom_y,
        u.mom_x * w.v,
        u.mom_y * w.v + w.p,
        (u.energy + w.p) * w.v,
    };
}

double sound_speed(const Primitive2D& w, double gamma) {
    return std::sqrt(gamma * w.p / w.rho);
}

double max_signal_speed_x(const Conservative2D& u, double gamma) {
    const Primitive2D w = conservative_to_primitive(u, gamma);
    return std::abs(w.u) + sound_speed(w, gamma);
}

double max_signal_speed_y(const Conservative2D& u, double gamma) {
    const Primitive2D w = conservative_to_primitive(u, gamma);
    return std::abs(w.v) + sound_speed(w, gamma);
}

Conservative2D rusanov_flux_x(const Conservative2D& left, const Conservative2D& right, double gamma) {
    // The retained 2D baseline uses local Lax--Friedrichs (Rusanov) dissipation.
    const Conservative2D f_left = physical_flux_x(left, gamma);
    const Conservative2D f_right = physical_flux_x(right, gamma);
    const double alpha = std::max(max_signal_speed_x(left, gamma), max_signal_speed_x(right, gamma));
    return 0.5 * (f_left + f_right) - 0.5 * alpha * (right - left);
}

Conservative2D rusanov_flux_y(const Conservative2D& bottom, const Conservative2D& top, double gamma) {
    // The y-direction route uses the corresponding normal signal speed.
    const Conservative2D f_bottom = physical_flux_y(bottom, gamma);
    const Conservative2D f_top = physical_flux_y(top, gamma);
    const double alpha = std::max(max_signal_speed_y(bottom, gamma), max_signal_speed_y(top, gamma));
    return 0.5 * (f_bottom + f_top) - 0.5 * alpha * (top - bottom);
}

double dx(const Grid2D& grid) {
    return (grid.x_max - grid.x_min) / static_cast<double>(grid.nx);
}

double dy(const Grid2D& grid) {
    return (grid.y_max - grid.y_min) / static_cast<double>(grid.ny);
}

double x_center(const Grid2D& grid, int i) {
    return grid.x_min + (static_cast<double>(i) + 0.5) * dx(grid);
}

double y_center(const Grid2D& grid, int j) {
    return grid.y_min + (static_cast<double>(j) + 0.5) * dy(grid);
}

}  // namespace fp
