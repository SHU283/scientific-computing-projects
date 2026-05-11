#include "numerical/advance.h"
#include "numerical/riemann_hll.h"
#include <algorithm>

// Swap the data stored in two fields
static inline void swap_fields(Field& A, Field& B)      
{
    std::swap(A.rho,  B.rho);
    std::swap(A.rhou, B.rhou);
    std::swap(A.rhov, B.rhov);
    std::swap(A.E,    B.E);
}

// First-order Godunov update
void advance_1st(Field& U, double dt, double dx, double dy, double gamma)  
{                                       
    // New field, stores Un+1
    Field Un(U.nx, U.ny, U.ng);    

    const double dtdx = dt / dx;   
    const double dtdy = dt / dy;

    // Only update the internal cells
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {  
            {
                // x-direction fluxes at i±1/2
                // Right face: (i, j) | (i+1, j)
                int idL = U.index(i,   j);
                int idR = U.index(i+1, j);
                Prim WL = cons_to_prim(U.rho[idL], U.rhou[idL], U.rhov[idL], U.E[idL], gamma);    
                Prim WR = cons_to_prim(U.rho[idR], U.rhou[idR], U.rhov[idR], U.E[idR], gamma);  
                Flux FxR = hll_flux_x(WL, WR, gamma); 

                // Left face: (i-1, j) | (i, j)
                int idLL = U.index(i-1, j);
                int idRR = U.index(i,   j);
                Prim WLL = cons_to_prim(U.rho[idLL], U.rhou[idLL], U.rhov[idLL], U.E[idLL], gamma);
                Prim WRR = cons_to_prim(U.rho[idRR], U.rhou[idRR], U.rhov[idRR], U.E[idRR], gamma);
                Flux FxL = hll_flux_x(WLL, WRR, gamma);

                // y-direction fluxes at j±1/2
                int idB = U.index(i, j-1);
                int idT = U.index(i, j);
                Prim WB = cons_to_prim(U.rho[idB], U.rhou[idB], U.rhov[idB], U.E[idB], gamma);
                Prim WT = cons_to_prim(U.rho[idT], U.rhou[idT], U.rhov[idT], U.E[idT], gamma);
                Flux GyB = hll_flux_y(WB, WT, gamma); // bottom face between (j-1) and j

                int idBB = U.index(i, j);
                int idTT = U.index(i, j+1);
                Prim WBB = cons_to_prim(U.rho[idBB], U.rhou[idBB], U.rhov[idBB], U.E[idBB], gamma);
                Prim WTT = cons_to_prim(U.rho[idTT], U.rhou[idTT], U.rhov[idTT], U.E[idTT], gamma);
                Flux GyT = hll_flux_y(WBB, WTT, gamma); // top face between j and (j+1)

                // Update conserved variables
                int id = U.index(i, j);

                Un.rho[id]  = U.rho[id]  - dtdx * (FxR.rho  - FxL.rho ) - dtdy * (GyT.rho  - GyB.rho );
                Un.rhou[id] = U.rhou[id] - dtdx * (FxR.rhou - FxL.rhou) - dtdy * (GyT.rhou - GyB.rhou);
                Un.rhov[id] = U.rhov[id] - dtdx * (FxR.rhov - FxL.rhov) - dtdy * (GyT.rhov - GyB.rhov);
                Un.E[id]    = U.E[id]    - dtdx * (FxR.E    - FxL.E   ) - dtdy * (GyT.E    - GyB.E   );
            }
        }
    }

    // Ghost does not update
    swap_fields(U, Un);
}
