#pragma once

#include <vector>

namespace fp {

// Conservative variables for the 1D Euler equations: U = (rho, rho*u, E).
struct Conservative1D {
    double rho = 0.0;
    double mom = 0.0;
    double energy = 0.0;
};

// Primitive variables: W = (rho, u, p).
struct Primitive1D {
    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;
};

// Cell-centered solution returned by a 1D solver run.
struct SimulationResult1D {
    std::vector<double> x;
    std::vector<Conservative1D> cells;
    double final_time = 0.0;
    int steps = 0;
};

Conservative1D operator+(const Conservative1D& a, const Conservative1D& b);
Conservative1D operator-(const Conservative1D& a, const Conservative1D& b);
Conservative1D operator*(double scalar, const Conservative1D& value);
Conservative1D operator*(const Conservative1D& value, double scalar);
Conservative1D operator/(const Conservative1D& value, double scalar);
Conservative1D& operator+=(Conservative1D& a, const Conservative1D& b);
Conservative1D& operator-=(Conservative1D& a, const Conservative1D& b);

// Conversions and flux utilities used by all 1D finite-volume methods.
Primitive1D conservative_to_primitive(const Conservative1D& u, double gamma);
Conservative1D primitive_to_conservative(const Primitive1D& w, double gamma);
Conservative1D physical_flux(const Conservative1D& u, double gamma);
double sound_speed(const Primitive1D& w, double gamma);

// Interface-speed and numerical-flux helpers for shock-capturing updates.
double max_wave_speed(const Conservative1D& left, const Conservative1D& right, double gamma);
Conservative1D hll_flux(const Conservative1D& left, const Conservative1D& right, double gamma);

// TVD limiter shared by MUSCL-style reconstructions.
double minmod(double a, double b);

}  // namespace fp
