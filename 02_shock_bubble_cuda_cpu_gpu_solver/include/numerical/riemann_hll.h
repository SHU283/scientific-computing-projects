#pragma once
#include "core/field.h"

// primitive
struct Prim {
    double rho;
    double u;
    double v;
    double p;
};

// conservative
struct Cons {
    double rho;
    double rhou;
    double rhov;
    double E;
}; 

// flux
struct Flux {
    double rho;
    double rhou;
    double rhov;
    double E;
};

// Convert conserved variables to primitive variables
Prim cons_to_prim(double rho, double rhou, double rhov, double E, double gamma);

// Physical flux in the x and y direction
Flux flux_x(const Prim& W, double gamma);       
Flux flux_y(const Prim& W, double gamma);       

// HLL numerical flux: x face (using u±c), y face (using v±c)
Flux hll_flux_x(const Prim& WL, const Prim& WR, double gamma);      
Flux hll_flux_y(const Prim& WL, const Prim& WR, double gamma);     
