#include "numerical/riemann_hll.h"
#include "physics/state_convert.h"
#include <algorithm>
#include <cmath>

// Cons to Prim variables
Prim cons_to_prim(double rho, double rhou, double rhov, double E, double gamma)     
{
    Prim W{};
    W.rho = rho;
    W.u   = rhou / rho;     
    W.v   = rhov / rho;
    W.p   = E_to_p(rho, rhou, rhov, E, gamma);
    return W;
}

// Physical flux in the x direction
Flux flux_x(const Prim& W, double /*gamma*/)       
{                                                   
    // F(U) in x-direction
    // [ rho*u,
    //   rho*u^2 + p,
    //   rho*u*v,
    //   (E+p)*u ]
    Flux F{};
    const double rho = W.rho;
    const double u   = W.u;
    const double v   = W.v;
    const double p   = W.p;

    F.rho  = rho * u;
    F.rhou = rho * u * u + p;
    F.rhov = rho * u * v;
    F.E    = 0.0;       // Fill in hll_flux
    return F;
}

// Physical flux in the y direction
Flux flux_y(const Prim& W, double /*gamma*/)        
{
    // G(U) in y-direction
    Flux G{};
    const double rho = W.rho;
    const double u   = W.u;
    const double v   = W.v;
    const double p   = W.p;

    G.rho  = rho * v;
    G.rhou = rho * u * v;
    G.rhov = rho * v * v + p;
    G.E    = 0.0;       // Fill in hll_flux
    return G;
}

// Prim to Cons 
static inline Cons prim_to_cons(const Prim& W, double gamma)        
{
    Cons U{};
    U.rho  = W.rho;
    U.rhou = W.rho * W.u;
    U.rhov = W.rho * W.v;
    U.E    = p_to_E(W.rho, W.u, W.v, W.p, gamma);
    return U;
}

// Cons to flux in x direction
static inline Flux cons_flux_x_full(const Cons& U, const Prim& W)       
{
    Flux F{};
    const double rho  = U.rho;
    const double rhou = U.rhou;
    const double rhov = U.rhov;
    const double E    = U.E; 
    const double u    = W.u;
    const double v    = W.v;
    const double p    = W.p;

    F.rho  = rhou;
    F.rhou = rhou * u + p;     // rho*u^2 + p
    F.rhov = rhov * u;         // rho*u*v
    F.E    = (E + p) * u;
    return F;
}

// Cons to flux in y direction
static inline Flux cons_flux_y_full(const Cons& U, const Prim& W)     
{
    Flux G{};
    const double rhou = U.rhou;
    const double rhov = U.rhov;
    const double E    = U.E;
    const double u    = W.u;
    const double v    = W.v;
    const double p    = W.p;

    G.rho  = rhov;
    G.rhou = rhou * v;         // rho*u*v
    G.rhov = rhov * v + p;     // rho*v^2 + p
    G.E    = (E + p) * v;
    return G;
}

// Compute the HLL flux
static inline Flux hll(const Cons& UL, const Cons& UR,    
                       const Flux& FL, const Flux& FR,
                       double SL, double SR)
{
    // Standard HLL flux:
    // if SL >= 0 -> FL
    // if SR <= 0 -> FR
    // else -> (SR*FL - SL*FR + SL*SR*(UR-UL)) / (SR-SL)
    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    const double inv = 1.0 / (SR - SL);

    Flux FH{};
    FH.rho  = (SR*FL.rho  - SL*FR.rho  + SL*SR*(UR.rho  - UL.rho )) * inv;
    FH.rhou = (SR*FL.rhou - SL*FR.rhou + SL*SR*(UR.rhou - UL.rhou)) * inv;
    FH.rhov = (SR*FL.rhov - SL*FR.rhov + SL*SR*(UR.rhov - UL.rhov)) * inv;
    FH.E    = (SR*FL.E    - SL*FR.E    + SL*SR*(UR.E    - UL.E   )) * inv;
    return FH;
}

// Calculate HLL flux in x direction
Flux hll_flux_x(const Prim& WL, const Prim& WR, double gamma)    
{
    // wave speeds in x: u ± c
    const double cL = std::sqrt(std::max(0.0, gamma * WL.p / WL.rho));
    const double cR = std::sqrt(std::max(0.0, gamma * WR.p / WR.rho));

    const double SL = std::min(WL.u - cL, WR.u - cR);
    const double SR = std::max(WL.u + cL, WR.u + cR);

    const Cons UL = prim_to_cons(WL, gamma);
    const Cons UR = prim_to_cons(WR, gamma);

    const Flux FL = cons_flux_x_full(UL, WL);
    const Flux FR = cons_flux_x_full(UR, WR);

    return hll(UL, UR, FL, FR, SL, SR);
}

// Calculate HLL flux in y direction
Flux hll_flux_y(const Prim& WL, const Prim& WR, double gamma)    
{
    // wave speeds in y: v ± c
    const double cL = std::sqrt(std::max(0.0, gamma * WL.p / WL.rho));
    const double cR = std::sqrt(std::max(0.0, gamma * WR.p / WR.rho));

    const double SL = std::min(WL.v - cL, WR.v - cR);
    const double SR = std::max(WL.v + cL, WR.v + cR);

    const Cons UL = prim_to_cons(WL, gamma);
    const Cons UR = prim_to_cons(WR, gamma);

    const Flux GL = cons_flux_y_full(UL, WL);
    const Flux GR = cons_flux_y_full(UR, WR);

    return hll(UL, UR, GL, GR, SL, SR);
}
