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
    const Conservative2D f_left = physical_flux_x(left, gamma);
    const Conservative2D f_right = physical_flux_x(right, gamma);
    const double alpha = std::max(max_signal_speed_x(left, gamma), max_signal_speed_x(right, gamma));
    return 0.5 * (f_left + f_right) - 0.5 * alpha * (right - left);
}

Conservative2D rusanov_flux_y(const Conservative2D& bottom, const Conservative2D& top, double gamma) {
    const Conservative2D f_bottom = physical_flux_y(bottom, gamma);
    const Conservative2D f_top = physical_flux_y(top, gamma);
    const double alpha = std::max(max_signal_speed_y(bottom, gamma), max_signal_speed_y(top, gamma));
    return 0.5 * (f_bottom + f_top) - 0.5 * alpha * (top - bottom);
}

Conservative2D hllc_flux_x(const Conservative2D& left, const Conservative2D& right, double gamma) {
    const Primitive2D wl = conservative_to_primitive(left, gamma);
    const Primitive2D wr = conservative_to_primitive(right, gamma);
    const double cl = sound_speed(wl, gamma);
    const double cr = sound_speed(wr, gamma);
    const double s_l = std::min(wl.u - cl, wr.u - cr);
    const double s_r = std::max(wl.u + cl, wr.u + cr);
    const Conservative2D f_l = physical_flux_x(left, gamma);
    const Conservative2D f_r = physical_flux_x(right, gamma);

    if (s_l >= 0.0) {
        return f_l;
    }
    if (s_r <= 0.0) {
        return f_r;
    }

    const double denom =
        wl.rho * (s_l - wl.u) - wr.rho * (s_r - wr.u);
    // Avoid amplifying roundoff when the contact-speed denominator is singular.
    if (std::abs(denom) < 1.0e-14) {
        return rusanov_flux_x(left, right, gamma);
    }

    const double s_m =
        (wr.p - wl.p + wl.rho * wl.u * (s_l - wl.u) -
         wr.rho * wr.u * (s_r - wr.u)) /
        denom;

    // The x-normal star state carries the original tangential y velocity.
    auto star_state = [&](const Conservative2D& state,
                          const Primitive2D& primitive,
                          double wave_speed) {
        const double rho_star =
            primitive.rho * (wave_speed - primitive.u) /
            (wave_speed - s_m);
        const double energy_per_mass = state.energy / primitive.rho;
        const double energy_star =
            rho_star *
            (energy_per_mass +
             (s_m - primitive.u) *
                 (s_m + primitive.p /
                            (primitive.rho * (wave_speed - primitive.u))));
        return Conservative2D{
            rho_star,
            rho_star * s_m,
            rho_star * primitive.v,
            energy_star,
        };
    };

    if (s_m >= 0.0) {
        const Conservative2D u_star_l = star_state(left, wl, s_l);
        return f_l + s_l * (u_star_l - left);
    }

    const Conservative2D u_star_r = star_state(right, wr, s_r);
    return f_r + s_r * (u_star_r - right);
}

Conservative2D hllc_flux_y(const Conservative2D& bottom, const Conservative2D& top, double gamma) {
    const Primitive2D wb = conservative_to_primitive(bottom, gamma);
    const Primitive2D wt = conservative_to_primitive(top, gamma);
    const double cb = sound_speed(wb, gamma);
    const double ct = sound_speed(wt, gamma);
    const double s_b = std::min(wb.v - cb, wt.v - ct);
    const double s_t = std::max(wb.v + cb, wt.v + ct);
    const Conservative2D f_b = physical_flux_y(bottom, gamma);
    const Conservative2D f_t = physical_flux_y(top, gamma);

    if (s_b >= 0.0) {
        return f_b;
    }
    if (s_t <= 0.0) {
        return f_t;
    }

    const double denom =
        wb.rho * (s_b - wb.v) - wt.rho * (s_t - wt.v);
    if (std::abs(denom) < 1.0e-14) {
        return rusanov_flux_y(bottom, top, gamma);
    }

    const double s_m =
        (wt.p - wb.p + wb.rho * wb.v * (s_b - wb.v) -
         wt.rho * wt.v * (s_t - wt.v)) /
        denom;

    // The y-normal construction mirrors the x flux and preserves tangential x velocity.
    auto star_state = [&](const Conservative2D& state,
                          const Primitive2D& primitive,
                          double wave_speed) {
        const double rho_star =
            primitive.rho * (wave_speed - primitive.v) /
            (wave_speed - s_m);
        const double energy_per_mass = state.energy / primitive.rho;
        const double energy_star =
            rho_star *
            (energy_per_mass +
             (s_m - primitive.v) *
                 (s_m + primitive.p /
                            (primitive.rho * (wave_speed - primitive.v))));
        return Conservative2D{
            rho_star,
            rho_star * primitive.u,
            rho_star * s_m,
            energy_star,
        };
    };

    if (s_m >= 0.0) {
        const Conservative2D u_star_b = star_state(bottom, wb, s_b);
        return f_b + s_b * (u_star_b - bottom);
    }

    const Conservative2D u_star_t = star_state(top, wt, s_t);
    return f_t + s_t * (u_star_t - top);
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
