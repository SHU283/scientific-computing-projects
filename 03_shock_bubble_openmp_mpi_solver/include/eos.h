#pragma once

// Conserved to primitive, calculate p
inline double pressure_from_cons(double rho, double rhou, double rhov, double E, double gamma) {        
    double u = rhou / rho;
    double v = rhov / rho;
    double kinetic = 0.5 * rho * (u*u + v*v);
    return (gamma - 1.0) * (E - kinetic);
}

// Primitive to conserved, calculate E
inline double energy_from_prim(double rho, double u, double v, double p, double gamma) {                
    double kinetic = 0.5 * rho * (u*u + v*v);
    double internal = p / (gamma - 1.0);
    return internal + kinetic;
}
