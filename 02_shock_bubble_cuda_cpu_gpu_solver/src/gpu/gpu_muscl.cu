#include <cuda_runtime.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "core/field.h"
#include "numerical/gpu_advance.h"
#include "core/gpu_field.h"

// Check
static void ck(cudaError_t e, const char* msg) {
    if (e != cudaSuccess) {
        std::cerr << "CUDA error: " << msg << " : "
                  << cudaGetErrorString(e) << "\n";
        std::exit(1);
    }
}

// Count the total number of cells on the host side
static inline int Ntot_host(const GPUField& U) { return U.pitch * U.height; }

// device helpers: Primitive / Conservative
struct Prim_d { double rho,u,v,p; };
struct Cons_d { double rho,rhou,rhov,E; };

// Compute pressure from cons
__device__ __forceinline__
double E_to_p_d(double rho, double rhou, double rhov, double E, double g){
    double u = rhou/rho, v = rhov/rho;
    return (g-1.0)*(E - 0.5*rho*(u*u + v*v));
}

// Compute total energy from prim
__device__ __forceinline__
double p_to_E_d(double rho, double u, double v, double p, double g){
    return p/(g-1.0) + 0.5*rho*(u*u + v*v);
}

// Cons to Prim
__device__ __forceinline__
Prim_d cons_to_prim_d(double rho, double rhou, double rhov, double E, double g){
    Prim_d W{};
    W.rho=rho;
    W.u=rhou/rho;
    W.v=rhov/rho;
    W.p=E_to_p_d(rho,rhou,rhov,E,g);
    return W;
}

// Prim to Cons
__device__ __forceinline__
Cons_d prim_to_cons_d(const Prim_d& q, double g){
    Cons_d U{};
    U.rho  = q.rho;
    U.rhou = q.rho*q.u;
    U.rhov = q.rho*q.v;
    U.E    = p_to_E_d(q.rho, q.u, q.v, q.p, g);
    return U;
}

// Calculate flux by Prim in x direction
__device__ __forceinline__
Cons_d flux_x_d(const Prim_d& q, double g){
    double E = p_to_E_d(q.rho, q.u, q.v, q.p, g);
    Cons_d F{};
    F.rho  = q.rho*q.u;
    F.rhou = q.rho*q.u*q.u + q.p;
    F.rhov = q.rho*q.u*q.v;
    F.E    = (E + q.p)*q.u;
    return F;
}

// Calculate flux by Prim in y direction
__device__ __forceinline__
Cons_d flux_y_d(const Prim_d& q, double g){
    double E = p_to_E_d(q.rho, q.u, q.v, q.p, g);
    Cons_d G{};
    G.rho  = q.rho*q.v;
    G.rhou = q.rho*q.u*q.v;
    G.rhov = q.rho*q.v*q.v + q.p;
    G.E    = (E + q.p)*q.v;
    return G;
}

// The speed of sound
__device__ __forceinline__
double sound_speed_d(double rho, double p, double g){
    return sqrt(fmax(0.0, g*p/rho));
}

// Perform operations on a 4-component Cons vector
__device__ __forceinline__
Cons_d addC(const Cons_d& A, const Cons_d& B){ return {A.rho+B.rho, A.rhou+B.rhou, A.rhov+B.rhov, A.E+B.E}; }
__device__ __forceinline__
Cons_d subC(const Cons_d& A, const Cons_d& B){ return {A.rho-B.rho, A.rhou-B.rhou, A.rhov-B.rhov, A.E-B.E}; }
__device__ __forceinline__
Cons_d mulC(const Cons_d& A, double s){ return {A.rho*s, A.rhou*s, A.rhov*s, A.E*s}; }

// Slope limiter
__device__ __forceinline__
double minmod_d(double a, double b){
    if (a*b <= 0.0) return 0.0;
    return (fabs(a) < fabs(b)) ? a : b;
}

// Perform operations on a 4-component Prim vector
__device__ __forceinline__
Prim_d addP(const Prim_d& A, const Prim_d& B){ return {A.rho+B.rho, A.u+B.u, A.v+B.v, A.p+B.p}; }
__device__ __forceinline__
Prim_d subP(const Prim_d& A, const Prim_d& B){ return {A.rho-B.rho, A.u-B.u, A.v-B.v, A.p-B.p}; }
__device__ __forceinline__
Prim_d mulP(const Prim_d& A, double s){ return {A.rho*s, A.u*s, A.v*s, A.p*s}; }

// Perform minmod operation on each of the four components of Prim variables
__device__ __forceinline__
Prim_d minmod_prim_d(const Prim_d& a, const Prim_d& b){
    return { minmod_d(a.rho,b.rho), minmod_d(a.u,b.u), minmod_d(a.v,b.v), minmod_d(a.p,b.p) };
}

// HLL flux in x direction (match muscl.cpp style: SL=min(u-c), SR=max(u+c))
__device__ __forceinline__
Cons_d hll_flux_x_d(const Prim_d& L, const Prim_d& R, double g){
    // Sound speeds and estimated wave speeds
    double cL = sound_speed_d(L.rho, L.p, g);
    double cR = sound_speed_d(R.rho, R.p, g);
    double SL = fmin(L.u - cL, R.u - cR);
    double SR = fmax(L.u + cL, R.u + cR);

    // left and right Prim to Cons, left and right physical fluxes
    Cons_d UL = prim_to_cons_d(L, g);
    Cons_d UR = prim_to_cons_d(R, g);
    Cons_d FL = flux_x_d(L, g);
    Cons_d FR = flux_x_d(R, g);

    // Upwind cases: use the left or right flux directly
    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    // Otherwise FHLL
    // (SR*FL - SL*FR + SL*SR*(UR-UL)) / (SR-SL)
    return mulC(
        addC(
            subC(mulC(FL, SR), mulC(FR, SL)),
            mulC(subC(UR, UL), SL*SR)
        ),
        1.0/(SR-SL)
    );
}

// HLL flux in y direction
__device__ __forceinline__
Cons_d hll_flux_y_d(const Prim_d& L, const Prim_d& R, double g){
    double cL = sound_speed_d(L.rho, L.p, g);
    double cR = sound_speed_d(R.rho, R.p, g);
    double SL = fmin(L.v - cL, R.v - cR);
    double SR = fmax(L.v + cL, R.v + cR);

    Cons_d UL = prim_to_cons_d(L, g);
    Cons_d UR = prim_to_cons_d(R, g);
    Cons_d GL = flux_y_d(L, g);
    Cons_d GR = flux_y_d(R, g);

    if (SL >= 0.0) return GL;
    if (SR <= 0.0) return GR;

    return mulC(
        addC(
            subC(mulC(GL, SR), mulC(GR, SL)),
            mulC(subC(UR, UL), SL*SR)
        ),
        1.0/(SR-SL)
    );
}

// temporary buffers for MUSCL (q_half + dq)
struct MusclTmp {
    double *qh_rho=nullptr, *qh_u=nullptr, *qh_v=nullptr, *qh_p=nullptr;
    double *dq_rho=nullptr, *dq_u=nullptr, *dq_v=nullptr, *dq_p=nullptr;
    int Nalloc=0;
};
static MusclTmp g_tmp;

// Ensure that the temporary array is large enough (perform malloc/free operations only on the host)
static void muscl_tmp_ensure(const GPUField& U){
    int N = Ntot_host(U);
    if (g_tmp.Nalloc == N) return;

    // Free current buffers
    auto fre = [&](double*& p){ if(p){ cudaFree(p); p=nullptr; } };
    fre(g_tmp.qh_rho); fre(g_tmp.qh_u); fre(g_tmp.qh_v); fre(g_tmp.qh_p);
    fre(g_tmp.dq_rho); fre(g_tmp.dq_u); fre(g_tmp.dq_v); fre(g_tmp.dq_p);

    // Allocate new buffers
    size_t bytes = sizeof(double)*(size_t)N;
    ck(cudaMalloc(&g_tmp.qh_rho, bytes), "malloc qh_rho");
    ck(cudaMalloc(&g_tmp.qh_u,   bytes), "malloc qh_u");
    ck(cudaMalloc(&g_tmp.qh_v,   bytes), "malloc qh_v");
    ck(cudaMalloc(&g_tmp.qh_p,   bytes), "malloc qh_p");

    ck(cudaMalloc(&g_tmp.dq_rho, bytes), "malloc dq_rho");
    ck(cudaMalloc(&g_tmp.dq_u,   bytes), "malloc dq_u");
    ck(cudaMalloc(&g_tmp.dq_v,   bytes), "malloc dq_v");
    ck(cudaMalloc(&g_tmp.dq_p,   bytes), "malloc dq_p");

    g_tmp.Nalloc = N;
}

// kernels: X sweep -------------------------------------------------------------------------

// 1) compute dq (minmod) and q_half (Hancock predictor) per cell in x direction
__global__ void x_predict_kernel(GPUField U, int cur,
                                 double lam, double g,
                                 double* qh_rho,double* qh_u,double* qh_v,double* qh_p,
                                 double* dq_rho,double* dq_u,double* dq_v,double* dq_p)
{
    // The thread is mapped to the grid point
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= U.pitch || j >= U.height) return;

    int ng=U.ng;
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;

    int id = i + U.pitch*j;

    // Not involved ghost / externally (with a zero slope, qh = current primitive)
    auto load_prim = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return cons_to_prim_d(U.rho[cur][k], U.rhou[cur][k], U.rhov[cur][k], U.E[cur][k], g);
    };

    // Set a default value for all the points first
    Prim_d q0 = load_prim(i,j);
    qh_rho[id]=q0.rho; qh_u[id]=q0.u; qh_v[id]=q0.v; qh_p[id]=q0.p;
    dq_rho[id]=0.0; dq_u[id]=0.0; dq_v[id]=0.0; dq_p[id]=0.0;

    // Exit immediately for ghost cells
    if (i<il || i>ir || j<jb || j>jt) return;

    // Primitive values of i-1, i, and i+1
    Prim_d qm1 = load_prim(i-1, j);
    Prim_d q0c = q0;
    Prim_d qp1 = load_prim(i+1, j);

    // Calculate the restricted slope dq (minmod limiter)
    Prim_d dq = minmod_prim_d(subP(q0c,qm1), subP(qp1,q0c));

    // Use dq to calculate the left and right interface values within the cell (still at tn)
    Prim_d q0L = subP(q0c, mulP(dq,0.5));
    Prim_d q0R = addP(q0c, mulP(dq,0.5));

    // Hancock predictor: First calculate the fluxes on the left and right sides 
    Cons_d F0L = flux_x_d(q0L, g);
    Cons_d F0R = flux_x_d(q0R, g);

    // Advance the cell-average to half time step
    Cons_d Uq0 = prim_to_cons_d(q0c, g);
    Cons_d Uq_half = subC(Uq0, mulC(subC(F0R,F0L), 0.5*lam));  // U - 0.5*lam*(F_R-F_L)

    // Cons to Prim, Save it in qh
    Prim_d qh = cons_to_prim_d(Uq_half.rho, Uq_half.rhou, Uq_half.rhov, Uq_half.E, g);

    qh_rho[id]=qh.rho; qh_u[id]=qh.u; qh_v[id]=qh.v; qh_p[id]=qh.p;
    dq_rho[id]=dq.rho; dq_u[id]=dq.u; dq_v[id]=dq.v; dq_p[id]=dq.p;
}

// 2) compute fluxes on-the-fly and update conserved into nxt
__global__ void x_update_kernel(GPUField U, int cur, int nxt,
                                double lam, double g,
                                const double* qh_rho,const double* qh_u,const double* qh_v,const double* qh_p,
                                const double* dq_rho,const double* dq_u,const double* dq_v,const double* dq_p)
{
    // Thread index and valid range check
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= U.pitch || j >= U.height) return;

    // Define the internal physical area
    int ng=U.ng;
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;

    int id = i + U.pitch*j;

    // ghost: copy through
    if (i<il || i>ir || j<jb || j>jt){
        U.rho[nxt][id]  = U.rho[cur][id];
        U.rhou[nxt][id] = U.rhou[cur][id];
        U.rhov[nxt][id] = U.rhov[cur][id];
        U.E[nxt][id]    = U.E[cur][id];
        return;
    }

    // Extract qh and dq (primitive structures) from the array
    auto load_qh = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return {qh_rho[k], qh_u[k], qh_v[k], qh_p[k]};
    };
    auto load_dq = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return {dq_rho[k], dq_u[k], dq_v[k], dq_p[k]};
    };

    // left interface flux F_{i-1/2}
    Cons_d Fm{};
    if (i == il){
        Prim_d q = load_qh(il, j);      // CPU: boundary interface uses first cell constant
        Fm = hll_flux_x_d(q, q, g);
    } else {
        Prim_d qLh = load_qh(i-1, j);
        Prim_d dqL = load_dq(i-1, j);
        Prim_d qRh = load_qh(i,   j);
        Prim_d dqR = load_dq(i,   j);

        Prim_d qL = addP(qLh, mulP(dqL,0.5)); // right side of left cell
        Prim_d qR = subP(qRh, mulP(dqR,0.5)); // left side of right cell
        Fm = hll_flux_x_d(qL, qR, g);
    }

    // right interface flux F_{i+1/2}
    Cons_d Fp{};
    if (i == ir){
        // Special handling of boundaries
        Prim_d q = load_qh(ir, j);
        Fp = hll_flux_x_d(q, q, g);
    } else {
        Prim_d qLh = load_qh(i,   j);
        Prim_d dqL = load_dq(i,   j);
        Prim_d qRh = load_qh(i+1, j);
        Prim_d dqR = load_dq(i+1, j);

        Prim_d qL = addP(qLh, mulP(dqL,0.5));
        Prim_d qR = subP(qRh, mulP(dqR,0.5));
        Fp = hll_flux_x_d(qL, qR, g);
    }

    // Update U_i: U^{n+1} = U^n - lam*(F_{i+1/2}-F_{i-1/2})
    double rho  = U.rho[cur][id];
    double rhou = U.rhou[cur][id];
    double rhov = U.rhov[cur][id];
    double E    = U.E[cur][id];

    U.rho[nxt][id]  = rho  - lam*(Fp.rho  - Fm.rho );
    U.rhou[nxt][id] = rhou - lam*(Fp.rhou - Fm.rhou);
    U.rhov[nxt][id] = rhov - lam*(Fp.rhov - Fm.rhov);
    U.E[nxt][id]    = E    - lam*(Fp.E    - Fm.E   );
}

// kernels: Y sweep -------------------------------------------------------------------------


// 1) compute dq (minmod) and q_half (Hancock predictor) per cell in y direction
__global__ void y_predict_kernel(GPUField U, int cur,
                                 double lam, double g,
                                 double* qh_rho,double* qh_u,double* qh_v,double* qh_p,
                                 double* dq_rho,double* dq_u,double* dq_v,double* dq_p)
{
    // The thread is mapped to the grid point
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= U.pitch || j >= U.height) return;

    int ng=U.ng;
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;

    int id = i + U.pitch*j;

    // Not involved ghost / externally (with a zero slope, qh = current primitive)
    auto load_prim = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return cons_to_prim_d(U.rho[cur][k], U.rhou[cur][k], U.rhov[cur][k], U.E[cur][k], g);
    };

    // Set a default value for all the points first
    Prim_d q0 = load_prim(i,j);
    qh_rho[id]=q0.rho; qh_u[id]=q0.u; qh_v[id]=q0.v; qh_p[id]=q0.p;
    dq_rho[id]=0.0; dq_u[id]=0.0; dq_v[id]=0.0; dq_p[id]=0.0;

    // Exit immediately for ghost cells
    if (i<il || i>ir || j<jb || j>jt) return;

    // Primitive values of j-1, j, and j+1 (y-direction stencil)
    Prim_d qm1 = load_prim(i, j-1);
    Prim_d q0c = q0;
    Prim_d qp1 = load_prim(i, j+1);

    // Calculate the restricted slope dq (minmod limiter)
    Prim_d dq = minmod_prim_d(subP(q0c,qm1), subP(qp1,q0c));

    // Use dq to calculate the bottom and top interface states within the cell (still at tn)
    Prim_d q0L = subP(q0c, mulP(dq,0.5));
    Prim_d q0R = addP(q0c, mulP(dq,0.5));

    // Hancock predictor: First calculate the fluxes on the left and right sides 
    Cons_d G0L = flux_y_d(q0L, g);
    Cons_d G0R = flux_y_d(q0R, g);

    // Advance the cell-average to half time step
    Cons_d Uq0 = prim_to_cons_d(q0c, g);
    Cons_d Uq_half = subC(Uq0, mulC(subC(G0R,G0L), 0.5*lam));

    // Cons to Prim, Save it in qh
    Prim_d qh = cons_to_prim_d(Uq_half.rho, Uq_half.rhou, Uq_half.rhov, Uq_half.E, g);

    qh_rho[id]=qh.rho; qh_u[id]=qh.u; qh_v[id]=qh.v; qh_p[id]=qh.p;
    dq_rho[id]=dq.rho; dq_u[id]=dq.u; dq_v[id]=dq.v; dq_p[id]=dq.p;
}

// 2) compute fluxes on-the-fly and update conserved into nxt
__global__ void y_update_kernel(GPUField U, int cur, int nxt,
                                double lam, double g,
                                const double* qh_rho,const double* qh_u,const double* qh_v,const double* qh_p,
                                const double* dq_rho,const double* dq_u,const double* dq_v,const double* dq_p)
{
    // Thread index and valid range check
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= U.pitch || j >= U.height) return;

    // Define the internal physical area
    int ng=U.ng;
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;

    int id = i + U.pitch*j;

    // ghost: copy through
    if (i<il || i>ir || j<jb || j>jt){
        U.rho[nxt][id]  = U.rho[cur][id];
        U.rhou[nxt][id] = U.rhou[cur][id];
        U.rhov[nxt][id] = U.rhov[cur][id];
        U.E[nxt][id]    = U.E[cur][id];
        return;
    }

    // Extract qh and dq (primitive structures) from the array
    auto load_qh = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return {qh_rho[k], qh_u[k], qh_v[k], qh_p[k]};
    };
    auto load_dq = [&](int ii,int jj)->Prim_d{
        int k = ii + U.pitch*jj;
        return {dq_rho[k], dq_u[k], dq_v[k], dq_p[k]};
    };

    // bottom interface G_{j-1/2}
    Cons_d Gm{};
    if (j == jb){
        Prim_d q = load_qh(i, jb);
        Gm = hll_flux_y_d(q, q, g);
    } else {
        Prim_d qBh = load_qh(i, j-1);
        Prim_d dqB = load_dq(i, j-1);
        Prim_d qTh = load_qh(i, j  );
        Prim_d dqT = load_dq(i, j  );

        Prim_d qB = addP(qBh, mulP(dqB,0.5));
        Prim_d qT = subP(qTh, mulP(dqT,0.5));
        Gm = hll_flux_y_d(qB, qT, g);
    }

    // top interface G_{j+1/2}
    Cons_d Gp{};
    if (j == jt){
        Prim_d q = load_qh(i, jt);
        Gp = hll_flux_y_d(q, q, g);
    } else {
        Prim_d qBh = load_qh(i, j  );
        Prim_d dqB = load_dq(i, j  );
        Prim_d qTh = load_qh(i, j+1);
        Prim_d dqT = load_dq(i, j+1);

        Prim_d qB = addP(qBh, mulP(dqB,0.5));
        Prim_d qT = subP(qTh, mulP(dqT,0.5));
        Gp = hll_flux_y_d(qB, qT, g);
    }

    // Update U_i: U^{n+1} = U^n - lam*(F_{j+1/2}-F_{j-1/2})
    double rho  = U.rho[cur][id];
    double rhou = U.rhou[cur][id];
    double rhov = U.rhov[cur][id];
    double E    = U.E[cur][id];

    U.rho[nxt][id]  = rho  - lam*(Gp.rho  - Gm.rho );
    U.rhou[nxt][id] = rhou - lam*(Gp.rhou - Gm.rhou);
    U.rhov[nxt][id] = rhov - lam*(Gp.rhov - Gm.rhov);
    U.E[nxt][id]    = E    - lam*(Gp.E    - Gm.E   );
}

// host wrappers: sweep_x / sweep_y / strang

// Perform a sweep in the x direction (with dt being either half-step or full-step)
static void sweep_x_gpu(GPUField& dU, double dt, double dx, double gamma){
    muscl_tmp_ensure(dU);

    // cur for reading, nxt for writing
    int cur = dU.cur;
    int nxt = 1 - dU.cur;

    double lam = dt / dx;

    // Configure the block/grid of the kernel
    dim3 block(16,16);
    dim3 grid((dU.pitch+block.x-1)/block.x,
              (dU.height+block.y-1)/block.y);

    // First, hancock predictor
    x_predict_kernel<<<grid,block>>>(dU, cur, lam, gamma,
        g_tmp.qh_rho,g_tmp.qh_u,g_tmp.qh_v,g_tmp.qh_p,
        g_tmp.dq_rho,g_tmp.dq_u,g_tmp.dq_v,g_tmp.dq_p);
    ck(cudaGetLastError(), "launch x_predict_kernel");

    // Calculate and update the interface HLL flux
    x_update_kernel<<<grid,block>>>(dU, cur, nxt, lam, gamma,
        g_tmp.qh_rho,g_tmp.qh_u,g_tmp.qh_v,g_tmp.qh_p,
        g_tmp.dq_rho,g_tmp.dq_u,g_tmp.dq_v,g_tmp.dq_p);
    ck(cudaGetLastError(), "launch x_update_kernel");

    // Synchronize and switch buffers
    ck(cudaDeviceSynchronize(), "sync sweep_x_gpu");
    dU.cur = nxt;
}

// Perform a sweep in the y direction (with dt being either half-step or full-step)
static void sweep_y_gpu(GPUField& dU, double dt, double dy, double gamma){
    muscl_tmp_ensure(dU);

    // cur for reading, nxt for writing
    int cur = dU.cur;
    int nxt = 1 - dU.cur;

    double lam = dt / dy;

    // Configure the block/grid of the kernel
    dim3 block(16,16);
    dim3 grid((dU.pitch+block.x-1)/block.x,
              (dU.height+block.y-1)/block.y);

    // First, hancock predictor
    y_predict_kernel<<<grid,block>>>(dU, cur, lam, gamma,
        g_tmp.qh_rho,g_tmp.qh_u,g_tmp.qh_v,g_tmp.qh_p,
        g_tmp.dq_rho,g_tmp.dq_u,g_tmp.dq_v,g_tmp.dq_p);
    ck(cudaGetLastError(), "launch y_predict_kernel");

    // Calculate and update the interface HLL flux
    y_update_kernel<<<grid,block>>>(dU, cur, nxt, lam, gamma,
        g_tmp.qh_rho,g_tmp.qh_u,g_tmp.qh_v,g_tmp.qh_p,
        g_tmp.dq_rho,g_tmp.dq_u,g_tmp.dq_v,g_tmp.dq_p);
    ck(cudaGetLastError(), "launch y_update_kernel");

    // Synchronize and switch buffers
    ck(cudaDeviceSynchronize(), "sync sweep_y_gpu");
    dU.cur = nxt;
}

void gpu_advance_2nd(GPUField& dU,
                                double dt, double dx, double dy, double gamma)
{
    // Use CUDA events to time the entire step
    cudaEvent_t s,e;
    cudaEventCreate(&s); cudaEventCreate(&e);
    cudaEventRecord(s);

    // Strang splitting
    sweep_x_gpu(dU, 0.5*dt, dx, gamma);
    sweep_y_gpu(dU, dt,     dy, gamma);
    sweep_x_gpu(dU, 0.5*dt, dx, gamma);

    // Record end time, synchronization, read ms
    cudaEventRecord(e);
    cudaEventSynchronize(e);

    float ms=0.f;
    cudaEventElapsedTime(&ms, s, e);
    cudaEventDestroy(s); cudaEventDestroy(e);

    // Update timing and step-count statistics
    gpu_timer_add_step_ms(ms);   
    gpu_timer_inc_steps(1);    
}

