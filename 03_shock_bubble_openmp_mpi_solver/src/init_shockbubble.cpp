#include "init.h"
#include "eos.h"
#include <cmath>
#include <cstdlib>

// Compute the post-shock state from the incoming Mach number
// 1: pre-shock, 2: post-shock
static void shock_post_state(double& rho2, double& u2, double& p2,
                             double rho1, double p1, double gamma, double M1)
{
    // Classical shock relations in the x-direction
    // First compute rho2/rho1 and p2/p1, then determine u2 from mass conservation
    double g = gamma;
    double M1s = M1*M1;

    double rho_ratio = ((g+1.0)*M1s) / ((g-1.0)*M1s + 2.0);
    double p_ratio   = 1.0 + 2.0*g/(g+1.0)*(M1s - 1.0);

    rho2 = rho_ratio * rho1;
    p2   = p_ratio * p1;

    // Sound speed and flow velocity
    double a1 = std::sqrt(g*p1/rho1);
    double u1 = M1 * a1;

    // Conservation of mass: rho1*u1 = rho2*u2
    u2 = (rho1*u1)/rho2;
}

// Initialize the bubble
void init_shock_bubble(Field& U, const Params& P)
{
    double dx = P.Lx / P.nx;
    double dy = P.Ly / P.ny;

    // pre-shock (region x > x_shock)
    double rho1 = P.rho_air;
    double p1   = P.p0;
    double u1   = 0.0;
    double v1   = 0.0;

    // post-shock (region x < x_shock)
    double rho2, p2, u2;
    shock_post_state(rho2, u2, p2, rho1, p1, P.gamma, P.M_shock);

    for (int j = 0; j < U.height; ++j) {
        for (int i = 0; i < U.pitch; ++i) {
            // cell center coordinates (only meaningful for interior)
            double x = (i - U.ng + 0.5) * dx;
            double y = (j - U.ng + 0.5) * dy;

            // base state = air, then bubble overwrites density
            bool shocked = (x < P.x_shock);

            double rho = shocked ? rho2 : rho1;
            double u   = shocked ? u2   : u1;
            double v   = v1;
            double p   = shocked ? p2   : p1;

            // bubble duplicated in y-direction for weak scaling
            bool in_bubble = false;
            for (int k = 0; k < std::max(1, P.weakN); ++k) {
                double yb_k = P.yb + k * 0.089;   // base Ly per tile
                double r = std::sqrt((x - P.xb)*(x - P.xb) + (y - yb_k)*(y - yb_k));
                if (r < P.rb) { in_bubble = true; break; }
            }
            if (in_bubble) rho = P.rho_he;

            int id = U.idx(i, j);
            U.rho[id]  = rho;
            U.rhou[id] = rho * u;
            U.rhov[id] = rho * v;
            U.E[id]    = energy_from_prim(rho, u, v, p, P.gamma);
        }
    }
}

void init_shock_bubble_local(Field& U, const Params& P, int j0_global)
{
    double dx = P.Lx / P.nx;
    double dy = P.Ly / P.ny;

    // pre-shock (region x > x_shock)
    double rho1 = P.rho_air;
    double p1   = P.p0;
    double u1   = 0.0;
    double v1   = 0.0;

    // post-shock (region x < x_shock)
    double rho2, p2, u2;
    shock_post_state(rho2, u2, p2, rho1, p1, P.gamma, P.M_shock);

    
    // U.ny is LOCAL interior height (without ghosts)
    // local interior j runs: j = ng ... ng + U.ny - 1
    // corresponding global interior j is: jg = j0_global + (j - ng)
    for (int j = 0; j < U.height; ++j) {
        // global interior j index (0..P.ny-1) for this local j
        int jg = j0_global + (j - U.ng);

        for (int i = 0; i < U.pitch; ++i) {
            int ig = i - U.ng; // global interior i index (0..P.nx-1) for this local i

            double x = (ig + 0.5) * dx;
            double y = (jg + 0.5) * dy;

            bool shocked = (x < P.x_shock);

            double rho = shocked ? rho2 : rho1;
            double u   = shocked ? u2   : u1;
            double v   = v1;
            double p   = shocked ? p2   : p1;

            bool in_bubble = false;
            for (int k = 0; k < std::max(1, P.weakN); ++k) {
                double yb_k = P.yb + k * 0.089;
                double r = std::sqrt((x - P.xb)*(x - P.xb) + (y - yb_k)*(y - yb_k));
                if (r < P.rb) { in_bubble = true; break; }
            }
            if (in_bubble) rho = P.rho_he;

            int id = U.idx(i, j);
            U.rho[id]  = rho;
            U.rhou[id] = rho * u;
            U.rhov[id] = rho * v;
            U.E[id]    = energy_from_prim(rho, u, v, p, P.gamma);
        }
    }
}
