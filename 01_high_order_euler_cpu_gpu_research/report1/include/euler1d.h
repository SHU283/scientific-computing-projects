#pragma once

#include <vector>

namespace fp {

// Conservative variables for the 1D Euler equations.
struct Conservative1D {
    double rho = 0.0;
    double mom = 0.0;
    double energy = 0.0;
};

// Primitive variables for the 1D Euler equations.
struct Primitive1D {
    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;
};

// Result object for a 1D simulation.
struct SimulationResult1D {
    // Cell-center coordinates.
    std::vector<double> x;
    // Conservative state stored in each grid cell.
    std::vector<Conservative1D> cells;
    // Final time reached by the simulation.
    double final_time = 0.0;
    // Total number of time steps taken.
    int steps = 0;
};

// Allow Conservative1D states to be combined with vector-like arithmetic.
Conservative1D operator+(const Conservative1D& a, const Conservative1D& b);
Conservative1D operator-(const Conservative1D& a, const Conservative1D& b);
Conservative1D operator*(double scalar, const Conservative1D& value);
Conservative1D operator*(const Conservative1D& value, double scalar);
Conservative1D operator/(const Conservative1D& value, double scalar);
Conservative1D& operator+=(Conservative1D& a, const Conservative1D& b);
Conservative1D& operator-=(Conservative1D& a, const Conservative1D& b);

// Convert conservative variables U=(rho,rho*u,E) to primitive variables W=(rho,u,p).
Primitive1D conservative_to_primitive(const Conservative1D& u, double gamma);
// Convert primitive variables W=(rho,u,p) to conservative variables U=(rho,rho*u,E).
Conservative1D primitive_to_conservative(const Primitive1D& w, double gamma);
// Physical Euler flux F(U).
Conservative1D physical_flux(const Conservative1D& u, double gamma);
// Thermodynamic sound speed.
double sound_speed(const Primitive1D& w, double gamma);
// Maximum characteristic speed across a two-state interface.
double max_wave_speed(const Conservative1D& left, const Conservative1D& right, double gamma);
// HLL numerical flux.
Conservative1D hll_flux(const Conservative1D& left, const Conservative1D& right, double gamma);
// Minmod limiter used for MUSCL slope reconstruction.
double minmod(double a, double b);

}  // namespace fp

