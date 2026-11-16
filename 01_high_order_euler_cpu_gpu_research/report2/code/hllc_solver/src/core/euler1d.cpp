#include "euler1d.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fp {

Conservative1D operator+(const Conservative1D& a, const Conservative1D& b) {
    return {a.rho + b.rho, a.mom + b.mom, a.energy + b.energy};
}

Conservative1D operator-(const Conservative1D& a, const Conservative1D& b) {
    return {a.rho - b.rho, a.mom - b.mom, a.energy - b.energy};
}

Conservative1D operator*(double scalar, const Conservative1D& value) {
    return {scalar * value.rho, scalar * value.mom, scalar * value.energy};
}

Conservative1D operator*(const Conservative1D& value, double scalar) {
    return scalar * value;
}

Conservative1D operator/(const Conservative1D& value, double scalar) {
    return {value.rho / scalar, value.mom / scalar, value.energy / scalar};
}

Conservative1D& operator+=(Conservative1D& a, const Conservative1D& b) {
    a.rho += b.rho;
    a.mom += b.mom;
    a.energy += b.energy;
    return a;
}

Conservative1D& operator-=(Conservative1D& a, const Conservative1D& b) {
    a.rho -= b.rho;
    a.mom -= b.mom;
    a.energy -= b.energy;
    return a;
}

Primitive1D conservative_to_primitive(const Conservative1D& u, double gamma) {
    if (u.rho <= 0.0) {
        throw std::runtime_error("Density must stay positive.");
    }

    const double velocity = u.mom / u.rho;
    const double kinetic = 0.5 * u.rho * velocity * velocity;
    const double pressure = (gamma - 1.0) * (u.energy - kinetic);

    if (pressure <= 0.0) {
        throw std::runtime_error("Pressure must stay positive.");
    }

    return {u.rho, velocity, pressure};
}

Conservative1D primitive_to_conservative(const Primitive1D& w, double gamma) {
    const double momentum = w.rho * w.u;
    const double internal = w.p / (gamma - 1.0);
    const double kinetic = 0.5 * w.rho * w.u * w.u;
    return {w.rho, momentum, internal + kinetic};
}

Conservative1D physical_flux(const Conservative1D& u, double gamma) {
    const Primitive1D w = conservative_to_primitive(u, gamma);
    return {
        u.mom,
        u.mom * w.u + w.p,
        (u.energy + w.p) * w.u,
    };
}

double sound_speed(const Primitive1D& w, double gamma) {
    return std::sqrt(gamma * w.p / w.rho);
}

double max_wave_speed(const Conservative1D& left, const Conservative1D& right, double gamma) {
    const Primitive1D wl = conservative_to_primitive(left, gamma);
    const Primitive1D wr = conservative_to_primitive(right, gamma);
    const double sl = std::abs(wl.u) + sound_speed(wl, gamma);
    const double sr = std::abs(wr.u) + sound_speed(wr, gamma);
    return std::max(sl, sr);
}

Conservative1D hll_flux(const Conservative1D& left, const Conservative1D& right, double gamma) {
    const Primitive1D wl = conservative_to_primitive(left, gamma);
    const Primitive1D wr = conservative_to_primitive(right, gamma);
    const double cl = sound_speed(wl, gamma);
    const double cr = sound_speed(wr, gamma);

    const double s_l = std::min(wl.u - cl, wr.u - cr);
    const double s_r = std::max(wl.u + cl, wr.u + cr);

    const Conservative1D f_l = physical_flux(left, gamma);
    const Conservative1D f_r = physical_flux(right, gamma);

    if (s_l >= 0.0) {
        return f_l;
    }
    if (s_r <= 0.0) {
        return f_r;
    }

    return (s_r * f_l - s_l * f_r + s_l * s_r * (right - left)) / (s_r - s_l);
}

// Resolve the left acoustic wave, contact wave, and right acoustic wave.
Conservative1D hllc_flux(const Conservative1D& left, const Conservative1D& right, double gamma) {
    const Primitive1D wl = conservative_to_primitive(left, gamma);
    const Primitive1D wr = conservative_to_primitive(right, gamma);
    const double cl = sound_speed(wl, gamma);
    const double cr = sound_speed(wr, gamma);

    const double s_l = std::min(wl.u - cl, wr.u - cr);
    const double s_r = std::max(wl.u + cl, wr.u + cr);
    const Conservative1D f_l = physical_flux(left, gamma);
    const Conservative1D f_r = physical_flux(right, gamma);

    if (s_l >= 0.0) {
        return f_l;
    }
    if (s_r <= 0.0) {
        return f_r;
    }

    const double denom =
        wl.rho * (s_l - wl.u) - wr.rho * (s_r - wr.u);
    // A nearly singular contact-speed formula is safer with the two-wave HLL flux.
    if (std::abs(denom) < 1.0e-14) {
        return hll_flux(left, right, gamma);
    }

    const double s_m =
        (wr.p - wl.p + wl.rho * wl.u * (s_l - wl.u) -
         wr.rho * wr.u * (s_r - wr.u)) /
        denom;

    // Rankine--Hugoniot relations determine the conservative state behind each wave.
    auto star_state = [&](const Conservative1D& state,
                          const Primitive1D& primitive,
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
        return Conservative1D{rho_star, rho_star * s_m, energy_star};
    };

    if (s_m >= 0.0) {
        const Conservative1D u_star_l = star_state(left, wl, s_l);
        return f_l + s_l * (u_star_l - left);
    }

    const Conservative1D u_star_r = star_state(right, wr, s_r);
    return f_r + s_r * (u_star_r - right);
}

double minmod(double a, double b) {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return (std::abs(a) < std::abs(b)) ? a : b;
}

}  // namespace fp
