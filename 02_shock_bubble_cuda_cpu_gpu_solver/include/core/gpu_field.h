#pragma once
#include "core/field.h"     

// GPU version of the field storage
struct GPUField {
    int nx, ny, ng, pitch, height;
    // Double buffers: one for the current state and one for the updated state
    double *rho[2], *rhou[2], *rhov[2], *E[2];
    int cur;
};

// Total number of cells including ghost cells
inline int gpu_size_tot(const GPUField& U){ return U.pitch * U.height; }    

// Memory management
// Allocate each variable on the GPU
void gpu_alloc(GPUField& dU, const Field& hU); 
// Release GPU memory
void gpu_free(GPUField& dU);  
// Copy conserved variables from CPU to GPU
void gpu_copy_c2g(GPUField& dU, const Field& hU);  
// Copy the results calculated on the GPU back to the CPU      
void gpu_copy_g2c(Field& hU, const GPUField& dU);    

// Field operations
// Fill ghost cells on the GPU
void gpu_bc(GPUField& dU);   
// Calculate CFL dt
double gpu_compute_dt_cfl(const GPUField& dU,double dx,double dy,double gamma,double CFL);      

// timer
// Number of GPU processing steps
long long gpu_get_step_calls();      
// Total time for BC
double gpu_get_bc_ms(); 
// Total time spent on the main steps                     
double gpu_get_step_ms();         
// Single-step execution time          
void gpu_timer_add_step_ms(double ms);      
// Increment the GPU step counter
void gpu_timer_inc_steps(long long n=1);    
