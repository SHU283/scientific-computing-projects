#include "hrsc.h"
#include "eos.h"
#include "bc.h"
#include "mpi_decomp.h"

#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

// Primitive variables
struct Prim {
    double rho, u, v, p;
};

// Conserved variable
struct Cons {
    double rho, rhou, rhov, E;
};

// Select a stable slope
static inline double minmod(double a, double b) {
    if (a*b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

// Convert conserved variables to primitive variables
static inline Prim cons_to_prim(const Cons& U, double gamma) {
    Prim W{};
    W.rho = std::max(U.rho, 1e-14);
    W.u   = U.rhou / W.rho;
    W.v   = U.rhov / W.rho;
    double p = pressure_from_cons(W.rho, U.rhou, U.rhov, U.E, gamma);
    W.p   = std::max(p, 1e-14);
    return W;
}

// Convert primitive variables to conserved variables
static inline Cons prim_to_cons(const Prim& W, double gamma) {
    Cons U{};
    double rho = std::max(W.rho, 1e-14);
    double p   = std::max(W.p,   1e-14);
    U.rho  = rho;
    U.rhou = rho * W.u;
    U.rhov = rho * W.v;
    U.E    = energy_from_prim(rho, W.u, W.v, p, gamma);
    return U;
}

// Compute x-direction flux from conserved variables
static inline void flux_x(const Cons& U, double gamma, double F[4]) {
    Prim W = cons_to_prim(U, gamma);
    double p = W.p;
    F[0] = U.rhou;
    F[1] = U.rhou * W.u + p;
    F[2] = U.rhou * W.v;
    F[3] = (U.E + p) * W.u;
}

// Compute y-direction flux from conserved variables
static inline void flux_y(const Cons& U, double gamma, double G[4]) {
    Prim W = cons_to_prim(U, gamma);
    double p = W.p;
    G[0] = U.rhov;
    G[1] = U.rhov * W.u;
    G[2] = U.rhov * W.v + p;
    G[3] = (U.E + p) * W.v;
}

// The HLL numerical flux on the x-direction interface
static inline void hll_flux_x(const Cons& UL, const Cons& UR, double gamma, double Fhll[4]) {

    // HLL requires the use of speed, sound velocity, pressure...
    Prim WL = cons_to_prim(UL, gamma);
    Prim WR = cons_to_prim(UR, gamma);

    // Sound speed (one on each side)
    double aL = std::sqrt(gamma * WL.p / WL.rho);
    double aR = std::sqrt(gamma * WR.p / WR.rho);

    // Left and right HLL wave speeds sL and sR
    double sL = std::min(WL.u - aL, WR.u - aR);
    double sR = std::max(WL.u + aL, WR.u + aR);

    // Store the left and right HLL fluxes
    double FL[4], FR[4];
    flux_x(UL, gamma, FL);
    flux_x(UR, gamma, FR);

    // The flux when the two wave velocities are the same
    if (sL >= 0.0) {
        for (int k=0;k<4;++k) Fhll[k] = FL[k];
        return;
    }
    if (sR <= 0.0) {
        for (int k=0;k<4;++k) Fhll[k] = FR[k];
        return;
    }

    // HLL
    for (int k=0;k<4;++k) {
        double ULk = (k==0?UL.rho:(k==1?UL.rhou:(k==2?UL.rhov:UL.E)));
        double URk = (k==0?UR.rho:(k==1?UR.rhou:(k==2?UR.rhov:UR.E)));
        Fhll[k] = (sR*FL[k] - sL*FR[k] + sL*sR*(URk - ULk)) / (sR - sL);
    }
}

// The HLL numerical flux on the y-direction interface
static inline void hll_flux_y(const Cons& UL, const Cons& UR, double gamma, double Ghll[4]) {
    Prim WL = cons_to_prim(UL, gamma);
    Prim WR = cons_to_prim(UR, gamma);

    double aL = std::sqrt(gamma * WL.p / WL.rho);
    double aR = std::sqrt(gamma * WR.p / WR.rho);

    double sL = std::min(WL.v - aL, WR.v - aR);
    double sR = std::max(WL.v + aL, WR.v + aR);

    double GL[4], GR[4];
    flux_y(UL, gamma, GL);
    flux_y(UR, gamma, GR);

    if (sL >= 0.0) {
        for (int k=0;k<4;++k) Ghll[k] = GL[k];
        return;
    }
    if (sR <= 0.0) {
        for (int k=0;k<4;++k) Ghll[k] = GR[k];
        return;
    }

    for (int k=0;k<4;++k) {
        double ULk = (k==0?UL.rho:(k==1?UL.rhou:(k==2?UL.rhov:UL.E)));
        double URk = (k==0?UR.rho:(k==1?UR.rhou:(k==2?UR.rhov:UR.E)));
        Ghll[k] = (sR*GL[k] - sL*GR[k] + sL*sR*(URk - ULk)) / (sR - sL);
    }
}

// Read and write conserved variables
static inline Cons get_cons(const Field& U, int id) {
    return Cons{U.rho[id], U.rhou[id], U.rhov[id], U.E[id]};
}
static inline void set_cons(Field& U, int id, const Cons& C) {
    U.rho[id]  = C.rho;
    U.rhou[id] = C.rhou;
    U.rhov[id] = C.rhov;
    U.E[id]    = C.E;
}

// thread-private workspace
struct WorkX {
    std::vector<double> F0, F1, F2, F3;
    std::vector<Cons> UL_edge, UR_edge;
    void ensure(int pitch) {
        if ((int)F0.size() != pitch+1) { F0.resize(pitch+1); F1.resize(pitch+1); F2.resize(pitch+1); F3.resize(pitch+1); }
        if ((int)UL_edge.size() != pitch) { UL_edge.resize(pitch); UR_edge.resize(pitch); }
    }
};

// Perform an MUSCL–Hancock finite volume update in the x-direction
void sweep_x(Field& U, double dt, double dx, double gamma) {
    Field Un(U.nx, U.ny, U.ng);
    Un.pitch = U.pitch; Un.height = U.height;
    // no full copy: only ensure storage
    if ((int)Un.rho.size() != U.pitch * U.height) {
        Un.rho.resize(U.pitch * U.height);
        Un.rhou.resize(U.pitch * U.height);
        Un.rhov.resize(U.pitch * U.height);
        Un.E.resize(U.pitch * U.height);
    }
    
    // Only write interior cells; ghosts will be handled by BC later
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        // per-thread reusable buffers
        static thread_local WorkX wk;
        wk.ensure(U.pitch);

        // Reconstruction + limiter + Hancock half-step prediction: Generating interface state
        for (int i = 1; i < U.pitch-1; ++i) {
            // Left, Center, Right primitive variables
            Prim Wm = cons_to_prim(get_cons(U, U.idx(i-1,j)), gamma);
            Prim W  = cons_to_prim(get_cons(U, U.idx(i  ,j)), gamma);
            Prim Wp = cons_to_prim(get_cons(U, U.idx(i+1,j)), gamma);

            // slope selection
            Prim dW{};
            dW.rho = minmod(W.rho - Wm.rho, Wp.rho - W.rho);
            dW.u   = minmod(W.u   - Wm.u,   Wp.u   - W.u);
            dW.v   = minmod(W.v   - Wm.v,   Wp.v   - W.v);
            dW.p   = minmod(W.p   - Wm.p,   Wp.p   - W.p);

            // Reconstruction values of the left and right interfaces within the same cell
            Prim WL{ W.rho - 0.5*dW.rho, W.u - 0.5*dW.u, W.v - 0.5*dW.v, W.p - 0.5*dW.p };
            Prim WR{ W.rho + 0.5*dW.rho, W.u + 0.5*dW.u, W.v + 0.5*dW.v, W.p + 0.5*dW.p };

            // Converting it into cons for flux
            Cons ULc = prim_to_cons(WL, gamma);
            Cons URc = prim_to_cons(WR, gamma);

            // Left and right original fluxes
            double FL[4], FR[4];
            flux_x(ULc, gamma, FL);
            flux_x(URc, gamma, FR);

            // Hancock half step
            double coef = (dt * 0.5) / dx;

            ULc.rho  += coef * (FL[0] - FR[0]);
            ULc.rhou += coef * (FL[1] - FR[1]);
            ULc.rhov += coef * (FL[2] - FR[2]);
            ULc.E    += coef * (FL[3] - FR[3]);

            URc.rho  += coef * (FL[0] - FR[0]);
            URc.rhou += coef * (FL[1] - FR[1]);
            URc.rhov += coef * (FL[2] - FR[2]);
            URc.E    += coef * (FL[3] - FR[3]);

            // Save results
            wk.UL_edge[i] = ULc;
            wk.UR_edge[i] = URc;
        }

        // Interface flux assignment
        for (int i = U.ng-1; i <= U.ng + U.nx - 1; ++i) {
            Cons ULc = wk.UR_edge[i];
            Cons URc = wk.UL_edge[i+1];
            double Fh[4];
            hll_flux_x(ULc, URc, gamma, Fh);
            wk.F0[i+1] = Fh[0];
            wk.F1[i+1] = Fh[1];
            wk.F2[i+1] = Fh[2];
            wk.F3[i+1] = Fh[3];
        }

        // Update each interior cell using the flux difference
        double coef = dt / dx;
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            int id = U.idx(i,j);
            Cons C = get_cons(U, id);

            C.rho  -= coef * (wk.F0[i+1] - wk.F0[i]);
            C.rhou -= coef * (wk.F1[i+1] - wk.F1[i]);
            C.rhov -= coef * (wk.F2[i+1] - wk.F2[i]);
            C.E    -= coef * (wk.F3[i+1] - wk.F3[i]);

            set_cons(Un, id, C);
        }
    }

    U.rho.swap(Un.rho);
    U.rhou.swap(Un.rhou);
    U.rhov.swap(Un.rhov);
    U.E.swap(Un.E);
}

// Same idea as WorkX, but for the y-direction
struct WorkY {
    std::vector<double> G0, G1, G2, G3;
    std::vector<Cons> UL_edge, UR_edge;
    void ensure(int height) {
        if ((int)G0.size() != height+1) { G0.resize(height+1); G1.resize(height+1); G2.resize(height+1); G3.resize(height+1); }
        if ((int)UL_edge.size() != height) { UL_edge.resize(height); UR_edge.resize(height); }
    }
};

void sweep_y(Field& U, double dt, double dy, double gamma) {
    Field Un(U.nx, U.ny, U.ng);
    Un.pitch = U.pitch; Un.height = U.height;
    // no full copy: only ensure storage
    if ((int)Un.rho.size() != U.pitch * U.height) {
        Un.rho.resize(U.pitch * U.height);
        Un.rhou.resize(U.pitch * U.height);
        Un.rhov.resize(U.pitch * U.height);
        Un.E.resize(U.pitch * U.height);
    }

    // we only write interior cells; ghosts will be handled by BC later
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = U.ng; i < U.ng + U.nx; ++i) {

        static thread_local WorkY wk;
        wk.ensure(U.height);

        for (int j = 1; j < U.height-1; ++j) {
            Prim Wm = cons_to_prim(get_cons(U, U.idx(i,j-1)), gamma);
            Prim W  = cons_to_prim(get_cons(U, U.idx(i,j  )), gamma);
            Prim Wp = cons_to_prim(get_cons(U, U.idx(i,j+1)), gamma);

            Prim dW{};
            dW.rho = minmod(W.rho - Wm.rho, Wp.rho - W.rho);
            dW.u   = minmod(W.u   - Wm.u,   Wp.u   - W.u);
            dW.v   = minmod(W.v   - Wm.v,   Wp.v   - W.v);
            dW.p   = minmod(W.p   - Wm.p,   Wp.p   - W.p);

            Prim WL{ W.rho - 0.5*dW.rho, W.u - 0.5*dW.u, W.v - 0.5*dW.v, W.p - 0.5*dW.p };
            Prim WR{ W.rho + 0.5*dW.rho, W.u + 0.5*dW.u, W.v + 0.5*dW.v, W.p + 0.5*dW.p };

            Cons ULc = prim_to_cons(WL, gamma);
            Cons URc = prim_to_cons(WR, gamma);

            double GL[4], GR[4];
            flux_y(ULc, gamma, GL);
            flux_y(URc, gamma, GR);

            double coef = (dt * 0.5) / dy;

            ULc.rho  += coef * (GL[0] - GR[0]);
            ULc.rhou += coef * (GL[1] - GR[1]);
            ULc.rhov += coef * (GL[2] - GR[2]);
            ULc.E    += coef * (GL[3] - GR[3]);

            URc.rho  += coef * (GL[0] - GR[0]);
            URc.rhou += coef * (GL[1] - GR[1]);
            URc.rhov += coef * (GL[2] - GR[2]);
            URc.E    += coef * (GL[3] - GR[3]);

            wk.UL_edge[j] = ULc;
            wk.UR_edge[j] = URc;
        }

        for (int j = U.ng-1; j <= U.ng + U.ny - 1; ++j) {
            Cons ULc = wk.UR_edge[j];
            Cons URc = wk.UL_edge[j+1];
            double Gh[4];
            hll_flux_y(ULc, URc, gamma, Gh);
            wk.G0[j+1] = Gh[0];
            wk.G1[j+1] = Gh[1];
            wk.G2[j+1] = Gh[2];
            wk.G3[j+1] = Gh[3];
        }

        double coef = dt / dy;
        for (int j = U.ng; j < U.ng + U.ny; ++j) {
            int id = U.idx(i,j);
            Cons C = get_cons(U, id);

            C.rho  -= coef * (wk.G0[j+1] - wk.G0[j]);
            C.rhou -= coef * (wk.G1[j+1] - wk.G1[j]);
            C.rhov -= coef * (wk.G2[j+1] - wk.G2[j]);
            C.E    -= coef * (wk.G3[j+1] - wk.G3[j]);

            set_cons(Un, id, C);
        }
    }

    U.rho.swap(Un.rho);
    U.rhou.swap(Un.rhou);
    U.rhov.swap(Un.rhov);
    U.E.swap(Un.E);
}

// Advance one step (Strang dimension splitting)
void advance_muscl_hancock(Field& U, double dt, const Params& P)
{
    double dx = P.Lx / P.nx;
    double dy = P.Ly / P.ny;

    sweep_x(U, 0.5*dt, dx, P.gamma);
    apply_bc_transmissive(U);          // Refresh ghost

    sweep_y(U,     dt, dy, P.gamma);
    apply_bc_transmissive(U);          

    sweep_x(U, 0.5*dt, dx, P.gamma);
    apply_bc_transmissive(U);         
}


