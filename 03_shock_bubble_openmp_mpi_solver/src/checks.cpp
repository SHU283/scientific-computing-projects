#include "checks.h"
#include "eos.h"
#include <cmath>
#include <iostream>

// Check whether the solution remains physically valid
bool check_physical(const Field& U, const Params& P, int step, double t)
{
    // Only scan the interior cells, excluding ghosts
    double rho_min = 1e30;   
    double p_min   = 1e30;
    int bad_count = 0;

    // Only scan the internal part, do not scan the ghost.
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            int id = U.idx(i, j);
            double rho = U.rho[id];

            // Check rho first to avoid dividing by non-positive density
            if (!std::isfinite(rho) || rho <= 0.0) {
                // Print the top five bad points
                if (bad_count < 5) {
                    std::cerr << "[BAD] step=" << step << " t=" << t
                              << " at (i,j)=(" << i << "," << j << ") "
                              << "rho=" << rho << " (skip p)\n";
                }
                bad_count++;
                continue;
            }

            double p = pressure_from_cons(rho, U.rhou[id], U.rhov[id], U.E[id], P.gamma);

            // If rho is valid, check p
            if (!std::isfinite(p) || p <= 0.0) {
                if (bad_count < 5) {
                    std::cerr << "[BAD] step=" << step << " t=" << t
                              << " at (i,j)=(" << i << "," << j << ") "
                              << "rho=" << rho << " p=" << p << "\n";
                }
                bad_count++;
            }

            // Update the minimum value
            if (rho < rho_min) rho_min = rho;
            if (p   < p_min)   p_min   = p;
        }
    }

    // open when need
    /*
    std::cout << "[CHECK] step=" << step << " t=" << t
              << " rho_min=" << rho_min << " p_min=" << p_min
              << " bad=" << bad_count << "\n";
    */

    return bad_count == 0;
}
