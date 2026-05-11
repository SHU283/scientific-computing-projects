#include "dt.h"
#include "eos.h"
#include "mpi_decomp.h"

#include <cmath>
#include <algorithm>
#include <limits>

double compute_dt_cfl(const Field& U, const Params& P)
{
    // Grid spacing
    double dx = P.Lx / P.nx;
    double dy = P.Ly / P.ny;

    // Initialize the minimum dt
    double dt_min = 1e30;

    // Determine the scanning range on y
    int j0 = U.ng;
    int j1 = U.ng + U.ny;

    // Exclude the boundary-adjacent interior rows in y and use only the
    // central local region for the CFL estimate
    if (U.ny > 2*U.ng) {        
        j0 = U.ng + U.ng;
        j1 = U.ng + U.ny - U.ng;
    }

    // Calculate the local dt for each cell in turn
    for (int j = j0; j < j1; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            int id = U.idx(i,j);

            // Recover variables and skip cells with very small density
            double rho = U.rho[id];
            if (rho <= 1e-14) continue;

            double u = U.rhou[id] / rho;
            double v = U.rhov[id] / rho;
            double p = pressure_from_cons(rho, U.rhou[id], U.rhov[id], U.E[id], P.gamma);
            p = std::max(p, 1e-14);

            // Calculate the speed of sound and the characteristic velocity
            double a = std::sqrt(P.gamma * p / rho);
            double sx = std::abs(u) + a;
            double sy = std::abs(v) + a;

            // Local CFL time step
            double dt_cell = 1.0 / (sx/dx + sy/dy);
            dt_min = std::min(dt_min, dt_cell);
        }
    }
    return P.CFL * dt_min;
}
