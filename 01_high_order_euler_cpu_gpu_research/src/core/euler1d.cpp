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
    // Convert U = (rho, rho*u, E) to W = (rho, u, p). Positivity checks are kept
    // close to the conversion because invalid density/pressure usually indicates
    // a failed numerical update near shocks or under too-large CFL settings.
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
    // Rebuild conservative variables from primitive input. The total energy is
    // internal energy plus kinetic energy, which is the quantity advanced by the
    // finite-volume update.
    const double momentum = w.rho * w.u;
    const double internal = w.p / (gamma - 1.0);
    const double kinetic = 0.5 * w.rho * w.u * w.u;
    return {w.rho, momentum, internal + kinetic};
}

Conservative1D physical_flux(const Conservative1D& u, double gamma) {
    // Physical Euler flux F(U) = (rho*u, rho*u^2+p, (E+p)u).
    // Numerical methods call this when constructing Riemann or predictor fluxes.
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
    // Conservative CFL estimates only need the fastest acoustic signal touching
    // an interface: |u| + c on either side.
    const Primitive1D wl = conservative_to_primitive(left, gamma);
    const Primitive1D wr = conservative_to_primitive(right, gamma);
    const double sl = std::abs(wl.u) + sound_speed(wl, gamma);
    const double sr = std::abs(wr.u) + sound_speed(wr, gamma);
    return std::max(sl, sr);
}

Conservative1D hll_flux(const Conservative1D& left, const Conservative1D& right, double gamma) {
    // HLL approximates the Riemann fan by two wave speeds. It is more robust than
    // a central flux for shocks while staying simple enough for CPU/GPU comparison.
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

double minmod(double a, double b) {
    // Return zero when neighboring slopes disagree in sign; otherwise keep the
    // smaller-magnitude slope. This is the basic TVD limiter used by MUSCL.
    if (a * b <= 0.0) {
        return 0.0;
    }
    return (std::abs(a) < std::abs(b)) ? a : b;
}

}  // namespace fp
