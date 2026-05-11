#pragma once
#include <cmath>

// Calculate the total energy density E from (ρ, u, v, p)
inline double p_to_E(double rho, double u, double v, double p, double g){        
    return p/(g-1.0) + 0.5*rho*(u*u + v*v);
}
// Calculate the pressure p from (ρ, ρu, ρv, E)
inline double E_to_p(double rho, double rhou, double rhov, double E, double g){  
    double u = rhou/rho, v = rhov/rho;
    return (g-1.0)*(E - 0.5*rho*(u*u + v*v));
}

/*
rho  : density
u, v : velocity components
p    : pressure
g    : ratio of specific heats
E    : total energy density
*/