#include "core/field.h"
#include "physics/state_convert.h"
#include <cmath>
#include <algorithm>
#include <limits>

double cfl_dt(const Field& U, double dx, double dy, double gamma, double CFL)
{
    // Set the initial value to infinity
    double dt_min = std::numeric_limits<double>::infinity();

    // Loop over all interior cells
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            
            int id = U.index(i, j);

            // Extract the conserved variables of the current cell
            double rho  = U.rho[id];
            double rhou = U.rhou[id];
            double rhov = U.rhov[id];
            double E    = U.E[id];

            // Safety check
            if (rho <= 1e-12) continue;
            double p = E_to_p(rho, rhou, rhov, E, gamma);
            if (p <= 1e-12) continue;

            // Velocity (u, v) and sound speed (c)
            double u = rhou / rho;
            double v = rhov / rho;
            double c = std::sqrt(gamma * p / rho);

            // Local time-step limits in x and y
            double max_speed_x = std::abs(u) + c;
            double max_speed_y = std::abs(v) + c;

            double dt_x = dx / max_speed_x;
            double dt_y = dy / max_speed_y;

            // Choose the smaller one of dt
            double dt_cell = std::min(dt_x, dt_y);

            // Update the smallest value of dt_min for the entire field
            if (dt_cell < dt_min) {
                dt_min = dt_cell;
            }
        }
    }

    // Safety check: If the field fails, provide a fallback value
    if (!std::isfinite(dt_min)) return 1e-12;

    return CFL * dt_min;
}

