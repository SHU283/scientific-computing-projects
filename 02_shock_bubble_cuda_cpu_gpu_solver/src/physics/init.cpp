#include "core/field.h"
#include "physics/state_convert.h" 
#include <cmath>

// Initialize the physical field
void init(Field& U, double Lx, double Ly,       
          double gamma, double p0,              
          double rho_air, double rho_he,       
          double xb, double yb, double rb,   
          double xshock)                    
{
    // The physical size of the cell
    double dx = Lx / U.nx;
    double dy = Ly / U.ny;

    // Shock Mach number
    const double M1 = 1.22;

    // Pre-shock state
    const double rho1 = rho_air;
    const double p1   = p0;
    const double a1 = std::sqrt(gamma * p1 / rho1);     // sound speed
    const double S  = M1 * a1;                          // shock speed

    // R – H jump condition leads to post-shock state
    const double p2   = p1 * (1.0 + 2.0 * gamma / (gamma + 1.0) * (M1*M1 - 1.0));
    const double rho2 = rho1 * ((gamma + 1.0) * M1*M1) / ((gamma - 1.0) * M1*M1 + 2.0);

    // R - H jump condition leads to the wave velocity u2
    const double u2 = S * (1.0 - rho1 / rho2);

    for (int j=0; j<U.height; ++j) {
        for (int i=0; i<U.pitch; ++i) {

            // Physical coordinates (ghost cells are included here and overwritten later by BC)
            double x = (i - U.ng + 0.5) * dx;
            double y = (j - U.ng + 0.5) * dy - 0.5*Ly; // Set y = 0 at the center line

            // Default state: Background air
            double rho = rho_air;
            double u   = 0.0, v = 0.0, p = p0;

            // Change the left side of the shock wave to the post-shock wave state
            if (x < xshock) {
                rho = rho2;
                u   = u2;
                p   = p2;
            }

            // Define helium bubble with a radius of rb
            double r2 = (x-xb)*(x-xb) + (y-yb)*(y-yb);
            if (r2 <= rb*rb && x >= xshock) {
                rho = rho_he;
                u   = 0.0;
                v   = 0.0;
                p   = p0;
            }

            // Prim to Cons
            int id = U.index(i,j);
            U.rho[id]  = rho;
            U.rhou[id] = rho * u;
            U.rhov[id] = rho * v;
            U.E[id]    = p_to_E(rho, u, v, p, gamma);
        }
    }
}
