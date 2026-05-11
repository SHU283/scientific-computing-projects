#include "numerical/advance.h"
#include "physics/state_convert.h"
#include <algorithm>
#include <cmath>
#include <vector>

// Primitive and conserved variables
struct Prim { double rho, u, v, p; };
struct Cons { double rho, rhou, rhov, E; };

// Read conserved variables from cell (i, j)
static inline Cons get_cons_ij(const Field& U, int i, int j) {            
    int id = U.index(i, j);
    return Cons{U.rho[id], U.rhou[id], U.rhov[id], U.E[id]};
}

// Calculate the speed of sound
static inline double sound_speed(double rho, double p, double gamma) {      
    return std::sqrt(gamma * p / rho);
}

//Cons to Prim
static inline Prim cons_to_prim(const Cons& U, double gamma) {      
    Prim q;
    q.rho = U.rho;
    q.u   = U.rhou / U.rho;
    q.v   = U.rhov / U.rho;
    q.p   = E_to_p(U.rho, U.rhou, U.rhov, U.E, gamma);
    return q;
}

// Prim to Cons
static inline Cons prim_to_cons(const Prim& q, double gamma) {      
    Cons U;
    U.rho  = q.rho;
    U.rhou = q.rho * q.u;
    U.rhov = q.rho * q.v;
    U.E    = p_to_E(q.rho, q.u, q.v, q.p, gamma);
    return U;
}

// Convert primitive variables to x-flux
static inline Cons flux_x(const Prim& q, double gamma) {            
    const double E = p_to_E(q.rho, q.u, q.v, q.p, gamma);
    return Cons{
        q.rho * q.u,
        q.rho * q.u * q.u + q.p,
        q.rho * q.u * q.v,
        (E + q.p) * q.u
    };
}

// Convert primitive variables to y-flux
static inline Cons flux_y(const Prim& q, double gamma) {            
    const double E = p_to_E(q.rho, q.u, q.v, q.p, gamma);
    return Cons{
        q.rho * q.v,
        q.rho * q.u * q.v,
        q.rho * q.v * q.v + q.p,
        (E + q.p) * q.v
    };
}

// Minmod slope limiter
static inline double minmod(double a, double b) {                   
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

// Apply minmod limiter to each primitive variable
static inline Prim minmod_prim(const Prim& a, const Prim& b) {      
    return Prim{
        minmod(a.rho, b.rho),
        minmod(a.u,   b.u),
        minmod(a.v,   b.v),
        minmod(a.p,   b.p)
    };
}

// Prim adding or subtracting or multiplication
static inline Prim sub(const Prim& A, const Prim& B) {              
    return Prim{A.rho - B.rho, A.u - B.u, A.v - B.v, A.p - B.p};
}
static inline Prim add(const Prim& A, const Prim& B) {
    return Prim{A.rho + B.rho, A.u + B.u, A.v + B.v, A.p + B.p};
}
static inline Prim mul(const Prim& A, double s) {
    return Prim{A.rho*s, A.u*s, A.v*s, A.p*s};
}

// Cons adding or subtracting or multiplication
static inline Cons subC(const Cons& A, const Cons& B) {             
    return Cons{A.rho-B.rho, A.rhou-B.rhou, A.rhov-B.rhov, A.E-B.E};
}
static inline Cons addC(const Cons& A, const Cons& B) {
    return Cons{A.rho+B.rho, A.rhou+B.rhou, A.rhov+B.rhov, A.E+B.E};
}
static inline Cons mulC(const Cons& A, double s) {
    return Cons{A.rho*s, A.rhou*s, A.rhov*s, A.E*s};
}

// Calculate the HLL flux of cell faces in x direction
static inline Cons hll_flux_x(const Prim& L, const Prim& R, double gamma) {     
    // Left-right sound speed
    const double cL = sound_speed(L.rho, L.p, gamma);           
    const double cR = sound_speed(R.rho, R.p, gamma);

    // Left-right wave speed
    const double SL = std::min(L.u - cL, R.u - cR);             
    const double SR = std::max(L.u + cL, R.u + cR);

    // Compute physical fluxes on both sides
    const Cons UL = prim_to_cons(L, gamma);                     
    const Cons UR = prim_to_cons(R, gamma);

    const Cons FL = flux_x(L, gamma);                           
    const Cons FR = flux_x(R, gamma);

    // When both waves are traveling to the right (left), the flux is FL (FR)
    if (SL >= 0.0) return FL;                                   
    if (SR <= 0.0) return FR;                                  

    // HLL flux for the intermediate wave configuration
    return mulC(
        addC(
            subC(mulC(FL, SR), mulC(FR, SL)),
            mulC(subC(UR, UL), SL*SR)
        ),
        1.0/(SR - SL)
    );
}

// Calculate the HLL flux of cell faces in y direction
static inline Cons hll_flux_y(const Prim& L, const Prim& R, double gamma) {
    // Left-right sound speed
    const double cL = sound_speed(L.rho, L.p, gamma);
    const double cR = sound_speed(R.rho, R.p, gamma);

    // Left-right wave speed
    const double SL = std::min(L.v - cL, R.v - cR);
    const double SR = std::max(L.v + cL, R.v + cR);

    // Compute physical fluxes on both sides
    const Cons UL = prim_to_cons(L, gamma);
    const Cons UR = prim_to_cons(R, gamma);

    const Cons GL = flux_y(L, gamma);
    const Cons GR = flux_y(R, gamma);

    // When both waves are traveling to the right (left), the flux is FL (FR)
    if (SL >= 0.0) return GL;
    if (SR <= 0.0) return GR;

    // HLL flux for the intermediate wave configuration
    return mulC(
        addC(
            subC(mulC(GL, SR), mulC(GR, SL)),
            mulC(subC(UR, UL), SL*SR)
        ),
        1.0/(SR - SL)
    );
}


// 1D MUSCL-Hancock sweep in x direction
// updates U in-place (conserved) but only interior cells
static void sweep_x_muscl(Field& U, double dt, double dx, double gamma) {
    const int ng = U.ng;
    const int nx = U.nx;
    const int ny = U.ny;

    const double lam = dt / dx;

    // Avoid read-write conflicts, record updated values
    std::vector<double> rho_new (U.rho.size());                                 
    std::vector<double> rhou_new(U.rhou.size());
    std::vector<double> rhov_new(U.rhov.size());
    std::vector<double> E_new   (U.E.size());

    // Copy the current field, including ghost cells
    rho_new  = U.rho;                           
    rhou_new = U.rhou;
    rhov_new = U.rhov;
    E_new    = U.E;
    
    // Process each interior row
    for (int j = ng; j < ng + ny; ++j) {     

        // nx + 1 interfaces in a row
        std::vector<Cons> F(nx + 1);            
        
        // each interface ii, construct the left and right states, calculate the hll flux
        for (int ii = 0; ii <= nx; ++ii) {          
            // The cell on the left (right) side, iL (iR)
            int iL = ng + ii - 1;                       
            int iR = ng + ii;                   

            // The primitive states on the left and right of the interface
            Prim qL, qR;        

            // Copy boundary cell, then calculate MUSCL-Hancock
            if (ii == 0) {                                  
                qL = cons_to_prim(get_cons_ij(U, ng, j), gamma);    
                qR = qL;
            } else if (ii == nx) {                          
                qL = cons_to_prim(get_cons_ij(U, ng + nx - 1, j), gamma);
                qR = qL;
            } else {
                // Get three primitives around the left cell (iL)
                Prim qm1 = cons_to_prim(get_cons_ij(U, iL - 1, j), gamma);          
                Prim q0  = cons_to_prim(get_cons_ij(U, iL, j), gamma);            
                Prim qp1 = cons_to_prim(get_cons_ij(U, iL + 1, j), gamma);          

                // Get three primitives around the right cell (iR)
                Prim rm1 = cons_to_prim(get_cons_ij(U, iR - 1, j), gamma);           
                Prim r0  = cons_to_prim(get_cons_ij(U, iR, j), gamma);
                Prim rp1 = cons_to_prim(get_cons_ij(U, iR + 1, j), gamma);

                // Select the left and right slopes
                Prim dqL = minmod_prim(sub(q0, qm1), sub(qp1, q0));        
                Prim dqR = minmod_prim(sub(r0, rm1), sub(rp1, r0));

                // Cell central value -> Cell left/right boundary value (left)
                Prim q0L = sub(q0, mul(dqL, 0.5));                        
                Prim q0R = add(q0, mul(dqL, 0.5));
                // The physical fluxes at the left and right boundary of the cell
                Cons F0L = flux_x(q0L, gamma);                          
                Cons F0R = flux_x(q0R, gamma);

                // Cell central value -> Cell left/right boundary value (right) 
                Prim r0L = sub(r0, mul(dqR, 0.5));                          
                Prim r0R = add(r0, mul(dqR, 0.5));
                Cons G0L = flux_x(r0L, gamma);                       
                Cons G0R = flux_x(r0R, gamma);

                // Hancock half-step prediction
                // Transform to Cons -> Half-step update -> Transform back to Prim -> Record
                Cons Uq0 = prim_to_cons(q0, gamma);                         
                Cons Ur0 = prim_to_cons(r0, gamma);

                Cons Uq_half = subC(Uq0, mulC(subC(F0R, F0L), 0.5 * lam));
                Cons Ur_half = subC(Ur0, mulC(subC(G0R, G0L), 0.5 * lam));

                Prim q_half = cons_to_prim(Uq_half, gamma);      
                Prim r_half = cons_to_prim(Ur_half, gamma);

                qL = add(q_half, mul(dqL, 0.5));                      
                qR = sub(r_half, mul(dqR, 0.5));   
            }

            // Store the HLL flux at this interface
            F[ii] = hll_flux_x(qL, qR, gamma);    
        }

        // Update each cell in this row using the flux difference
        for (int ii = 0; ii < nx; ++ii) {                              
            int i = ng + ii;
            int id = U.index(i, j);

            // Store the old value and the new value
            Cons Uold{U.rho[id], U.rhou[id], U.rhov[id], U.E[id]};    
            Cons Unew = subC(Uold, mulC(subC(F[ii+1], F[ii]), lam));

            rho_new[id]  = Unew.rho;
            rhou_new[id] = Unew.rhou;
            rhov_new[id] = Unew.rhov;
            E_new[id]    = Unew.E;
        }
    }

    // Replace the old field with the updated field
    U.rho.swap(rho_new);                         
    U.rhou.swap(rhou_new);
    U.rhov.swap(rhov_new);
    U.E.swap(E_new);
}

// Repeat the same process in the y direction
static void sweep_y_muscl(Field& U, double dt, double dy, double gamma) {  
    const int ng = U.ng;
    const int nx = U.nx;
    const int ny = U.ny;

    const double lam = dt / dy;

    std::vector<double> rho_new (U.rho.size());
    std::vector<double> rhou_new(U.rhou.size());
    std::vector<double> rhov_new(U.rhov.size());
    std::vector<double> E_new   (U.E.size());

    rho_new  = U.rho;
    rhou_new = U.rhou;
    rhov_new = U.rhov;
    E_new    = U.E;

    for (int i = ng; i < ng + nx; ++i) {
        std::vector<Cons> G(ny + 1);

        for (int jj = 0; jj <= ny; ++jj) {
            int jB = ng + jj - 1;
            int jT = ng + jj;

            Prim qB, qT;

            if (jj == 0) {
                qB = cons_to_prim(get_cons_ij(U, i, ng), gamma);
                qT = qB;
            } else if (jj == ny) {
                qB = cons_to_prim(get_cons_ij(U, i, ng + ny - 1), gamma);
                qT = qB;
            } else {
                Prim qm1 = cons_to_prim(get_cons_ij(U, i, jB - 1), gamma);
                Prim q0  = cons_to_prim(get_cons_ij(U, i, jB), gamma);
                Prim qp1 = cons_to_prim(get_cons_ij(U, i, jB + 1), gamma);

                Prim rm1 = cons_to_prim(get_cons_ij(U, i, jT - 1), gamma);
                Prim r0  = cons_to_prim(get_cons_ij(U, i, jT), gamma);
                Prim rp1 = cons_to_prim(get_cons_ij(U, i, jT + 1), gamma);

                Prim dqB = minmod_prim(sub(q0, qm1), sub(qp1, q0));
                Prim dqT = minmod_prim(sub(r0, rm1), sub(rp1, r0));

                Prim q0L = sub(q0, mul(dqB, 0.5));
                Prim q0R = add(q0, mul(dqB, 0.5));
                Cons G0L = flux_y(q0L, gamma);
                Cons G0R = flux_y(q0R, gamma);

                Prim r0L = sub(r0, mul(dqT, 0.5));
                Prim r0R = add(r0, mul(dqT, 0.5));
                Cons H0L = flux_y(r0L, gamma);
                Cons H0R = flux_y(r0R, gamma);

                Cons Uq0 = prim_to_cons(q0, gamma);
                Cons Ur0 = prim_to_cons(r0, gamma);

                Cons Uq_half = subC(Uq0, mulC(subC(G0R, G0L), 0.5 * lam));
                Cons Ur_half = subC(Ur0, mulC(subC(H0R, H0L), 0.5 * lam));

                Prim q_half = cons_to_prim(Uq_half, gamma);
                Prim r_half = cons_to_prim(Ur_half, gamma);

                qB = add(q_half, mul(dqB, 0.5));
                qT = sub(r_half, mul(dqT, 0.5));
            }

            G[jj] = hll_flux_y(qB, qT, gamma);
        }

        for (int jj = 0; jj < ny; ++jj) {
            int j = ng + jj;
            int id = U.index(i, j);

            Cons Uold{U.rho[id], U.rhou[id], U.rhov[id], U.E[id]};
            Cons Unew = subC(Uold, mulC(subC(G[jj+1], G[jj]), lam));

            rho_new[id]  = Unew.rho;
            rhou_new[id] = Unew.rhou;
            rhov_new[id] = Unew.rhov;
            E_new[id]    = Unew.E; 
        }
    }

    U.rho.swap(rho_new);
    U.rhou.swap(rhou_new);
    U.rhov.swap(rhov_new);
    U.E.swap(E_new);
}

// Advance one step using Strang splitting
void advance_2nd(Field& U, double dt, double dx, double dy, double gamma) {     

    // x half-step
    sweep_x_muscl(U, 0.5 * dt, dx, gamma);

    // y full-step
    sweep_y_muscl(U, dt, dy, gamma);

    // x half-step
    sweep_x_muscl(U, 0.5 * dt, dx, gamma);
}
