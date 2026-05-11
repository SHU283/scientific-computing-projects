#include "core/field.h"      
#include <cuda_runtime.h>
#include <iostream>
#include <cmath>
#include <algorithm> 
#include "core/gpu_field.h"

// Global counter: cumulative time + call count
static double g_kernel_ms = 0.0;                                                      
static long long g_step_calls = 0;                                                 

// Increment the GPU step counter
void gpu_step_count_hook(){                                                             
    g_step_calls++;
}

// Check
static void ck(cudaError_t e, const char* msg) {                                        
    if (e != cudaSuccess) {
        std::cerr << "CUDA error: " << msg << " : "
                  << cudaGetErrorString(e) << "\n";
        std::exit(1);
    }
}

// Two types of variables: primitive variables vs conserved variables
struct Prim_d { double rho,u,v,p; };                                                   
struct Flux_d { double rho,rhou,rhov,E; };

// Calculate p using E
__device__ __forceinline__
double E_to_p_d(double rho, double rhou, double rhov, double E, double g){              
    double u = rhou/rho, v = rhov/rho;
    return (g-1.0)*(E - 0.5*rho*(u*u + v*v));
}

// Calculate E using p
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

// Cons to flux in x direction
__device__ __forceinline__
Flux_d cons_flux_x_full_d(double rho, double rhou, double rhov, double E,               
                          const Prim_d& W)
{
    Flux_d F{};
    double u=W.u, p=W.p;
    F.rho  = rhou;
    F.rhou = rhou*u + p;
    F.rhov = rhov*u;
    F.E    = (E + p)*u;
    return F;
}

// Cons to flux in y direction
__device__ __forceinline__
Flux_d cons_flux_y_full_d(double rho, double rhou, double rhov, double E,               
                          const Prim_d& W)
{
    Flux_d G{};
    double v=W.v, p=W.p;
    G.rho  = rhov;
    G.rhou = rhou*v;
    G.rhov = rhov*v + p;
    G.E    = (E + p)*v;
    return G;
}

// Core HLL flux formula
__device__ __forceinline__
Flux_d hll_core_d(double UL_rho, double UL_rhou, double UL_rhov, double UL_E,           
                  double UR_rho, double UR_rhou, double UR_rhov, double UR_E,
                  const Flux_d& FL, const Flux_d& FR,
                  double SL, double SR)
{
    // If both waves are moving to the right (left), then use the left (right) flux
    if (SL >= 0.0) return FL;                                                           
    if (SR <= 0.0) return FR;     
    // Otherwise, FHLL                                                     
    double inv = 1.0/(SR-SL);                                                         
    Flux_d FH{};
    FH.rho  = (SR*FL.rho  - SL*FR.rho  + SL*SR*(UR_rho  - UL_rho )) * inv;              
    FH.rhou = (SR*FL.rhou - SL*FR.rhou + SL*SR*(UR_rhou - UL_rhou)) * inv;
    FH.rhov = (SR*FL.rhov - SL*FR.rhov + SL*SR*(UR_rhov - UL_rhov)) * inv;
    FH.E    = (SR*FL.E    - SL*FR.E    + SL*SR*(UR_E    - UL_E   )) * inv;
    return FH;
}

// Calculate the HLL flux of the interface in the x direction
__device__ __forceinline__
Flux_d hll_flux_x_d(const Prim_d& WL, const Prim_d& WR, double g)                      
{
    // Calculate the speeds of sound cL and cR
    double cL = sqrt(fmax(0.0, g * WL.p / WL.rho));                                     
    double cR = sqrt(fmax(0.0, g * WR.p / WR.rho));
    // Select wave velocity SL, SR (conservative)
    double SL = fmin(WL.u - cL, WR.u - cR);                                
    double SR = fmax(WL.u + cL, WR.u + cR);

    // Prim to Cons for UL
    double UL_rho  = WL.rho;                                                          
    double UL_rhou = WL.rho * WL.u;
    double UL_rhov = WL.rho * WL.v;
    double UL_E    = p_to_E_d(WL.rho, WL.u, WL.v, WL.p, g);

    // Prim to Cons for UR
    double UR_rho  = WR.rho;                                                           
    double UR_rhou = WR.rho * WR.u;
    double UR_rhov = WR.rho * WR.v;
    double UR_E    = p_to_E_d(WR.rho, WR.u, WR.v, WR.p, g);

    // Left and right physical fluxes FL, FR
    Flux_d FL = cons_flux_x_full_d(UL_rho, UL_rhou, UL_rhov, UL_E, WL);               
    Flux_d FR = cons_flux_x_full_d(UR_rho, UR_rhou, UR_rhov, UR_E, WR);

    // Determine interface flux
    return hll_core_d(UL_rho,UL_rhou,UL_rhov,UL_E,                                    
                      UR_rho,UR_rhou,UR_rhov,UR_E,
                      FL,FR,SL,SR);
}

// Calculate the HLL flux of the interface in the y direction
__device__ __forceinline__
Flux_d hll_flux_y_d(const Prim_d& WL, const Prim_d& WR, double g)                    
{
    double cL = sqrt(fmax(0.0, g * WL.p / WL.rho));
    double cR = sqrt(fmax(0.0, g * WR.p / WR.rho));
    double SL = fmin(WL.v - cL, WR.v - cR);
    double SR = fmax(WL.v + cL, WR.v + cR);

    double UL_rho  = WL.rho;
    double UL_rhou = WL.rho * WL.u;
    double UL_rhov = WL.rho * WL.v;
    double UL_E    = p_to_E_d(WL.rho, WL.u, WL.v, WL.p, g);

    double UR_rho  = WR.rho;
    double UR_rhou = WR.rho * WR.u;
    double UR_rhov = WR.rho * WR.v;
    double UR_E    = p_to_E_d(WR.rho, WR.u, WR.v, WR.p, g);

    Flux_d GL = cons_flux_y_full_d(UL_rho, UL_rhou, UL_rhov, UL_E, WL);
    Flux_d GR = cons_flux_y_full_d(UR_rho, UR_rhou, UR_rhov, UR_E, WR);

    return hll_core_d(UL_rho,UL_rhou,UL_rhov,UL_E,
                      UR_rho,UR_rhou,UR_rhov,UR_E,
                      GL,GR,SL,SR);
}

// Update from Un to Un+1
__global__ void step_hll_1st_kernel(GPUField U,                                      
                                   int cur, int nxt,
                                   double dt, double dx, double dy, double g)
{
    // Each thread corresponds to one grid point
    int i = blockIdx.x * blockDim.x + threadIdx.x;                                    
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    // Threads beyond the specified range do not need to be considered
    if (i >= U.pitch || j >= U.height) return;                                      

    // Physical area scope
    int ng=U.ng;                                                                      
    int il=ng, ir=ng+U.nx-1;                                                      
    int jb=ng, jt=ng+U.ny-1;

    // Map (i, j) to a one-dimensional id
    int id = i + U.pitch*j;                                                     

    // Ghost area does not update: Simply copy Un → Un+1
    if (i<il || i>ir || j<jb || j>jt){                                           
        U.rho[nxt][id]  = U.rho[cur][id];
        U.rhou[nxt][id] = U.rhou[cur][id];
        U.rhov[nxt][id] = U.rhov[cur][id];
        U.E[nxt][id]    = U.E[cur][id];
        return;
    }

    // Read conserved variables from the current buffer and convert to primitive variables
    auto prim_at = [&](int ii, int jj) {                                              
        int k = ii + U.pitch * jj;                                                   
        return cons_to_prim_d(U.rho[cur][k], U.rhou[cur][k],
                              U.rhov[cur][k], U.E[cur][k], g);
    };

    // Left and right states in the x direction
    Prim_d WL  = prim_at(i,   j);                                                   
    Prim_d WR  = prim_at(i+1, j);                                                       
    Prim_d WLL = prim_at(i-1, j);                                                       

    // Left and right states in the y direction
    Prim_d WB  = prim_at(i, j-1);                                                     
    Prim_d WT  = prim_at(i, j);
    Prim_d WTT = prim_at(i, j+1);

    // Left and right interface flux FH
    Flux_d FxR = hll_flux_x_d(WL, WR, g);                                               
    Flux_d FxL = hll_flux_x_d(WLL, WL, g);                                           

    // Top and bottom interface flux FH
    Flux_d GyB = hll_flux_y_d(WB, WT, g);                                              
    Flux_d GyT = hll_flux_y_d(WT, WTT, g);                                        

    double dtdx = dt/dx;                                                       
    double dtdy = dt/dy;

    // Read cons (Un) of the current grid point
    double rho  = U.rho[cur][id];                                               
    double rhou = U.rhou[cur][id];
    double rhov = U.rhov[cur][id];
    double E    = U.E[cur][id];

    // Advance each variable to Un+1
    U.rho[nxt][id]  = rho  - dtdx*(FxR.rho  - FxL.rho ) - dtdy*(GyT.rho  - GyB.rho );
    U.rhou[nxt][id] = rhou - dtdx*(FxR.rhou - FxL.rhou) - dtdy*(GyT.rhou - GyB.rhou);
    U.rhov[nxt][id] = rhov - dtdx*(FxR.rhov - FxL.rhov) - dtdy*(GyT.rhov - GyB.rhov);
    U.E[nxt][id]    = E    - dtdx*(FxR.E    - FxL.E   ) - dtdy*(GyT.E    - GyB.E   );
}

// Advance one time step (1st)
void gpu_advance_1st(GPUField& dU,double dt, double dx, double dy, double gamma)        
{
    // cur for reading, nxt for writing
    int cur = dU.cur;                                                                   
    int nxt = 1 - dU.cur;                                                           

    // CUDA event: Timing on GPU side
    cudaEvent_t s,e;                                                               
    cudaEventCreate(&s); cudaEventCreate(&e);

    // Configure thread block and grid
    dim3 block(16,16);                                                                 
    dim3 grid((dU.pitch+block.x-1)/block.x,                                      
              (dU.height+block.y-1)/block.y);

    // Run the kernel and time it
    cudaEventRecord(s);                                                             
    step_hll_1st_kernel<<<grid,block>>>(dU, cur, nxt, dt, dx, dy, gamma);
    cudaEventRecord(e);

    // Error checking + Wait for GPU to finish running
    ck(cudaGetLastError(), "launch step_hll_1st_kernel");                          
    ck(cudaEventSynchronize(e), "sync step event");                            

    // Calculate the time consumed by the kernel and accumulate it into the global statistics
    float ms=0.f;                                                                 
    cudaEventElapsedTime(&ms, s, e);
    g_kernel_ms += ms;

    // Clear event resources
    cudaEventDestroy(s); cudaEventDestroy(e);                                       

    // Switch buffer + Count the number of steps
    dU.cur = nxt;                                                                      
    g_step_calls++;                                                                    
}

// Return the number of GPU update steps
long long gpu_get_step_calls(){ return g_step_calls; }                         

// Each unit sets an upper limit for dt according to the CFL
__device__ __forceinline__
double safe_dt_cell(double rho, double rhou, double rhov, double E,                    
                    double dx, double dy, double g)
{
    if (rho <= 0.0) return 1e30;

    // Cons to Prim
    double u = rhou / rho;                                                              
    double v = rhov / rho;
    double p = E_to_p_d(rho, rhou, rhov, E, g);
    if (p <= 0.0) return 1e30;

    // Sound speed c and the maximum characteristic velocity
    double c = sqrt(fmax(0.0, g * p / rho));                                       
    double smax_x = fabs(u) + c;                                                   
    double smax_y = fabs(v) + c;

    // Safety check
    if (smax_x <= 0.0 || smax_y <= 0.0) return 1e30;

    // CFL dt: Calculate dx/velocity and dy/velocity separately, and choose the smaller one
    double dt_x = dx / smax_x;
    double dt_y = dy / smax_y;
    return fmin(dt_x, dt_y);                                                     
}

// Cumulative kernel milliseconds + manual addition time + manual addition steps
double gpu_get_step_ms(){ return g_kernel_ms; }                                         
void gpu_timer_add_step_ms(double ms){ g_kernel_ms += ms; }                             
void gpu_timer_inc_steps(long long n){ g_step_calls += n; }                        
