#include "core/gpu_field.h"
#include <cuda_runtime.h>
#include <iostream>
#include <cstdlib>

// Global cumulative BC time
static double g_bc_ms = 0.0;            

// GPU error detection
static void ck(cudaError_t e, const char* msg){             
    if(e!=cudaSuccess){
        std::cerr<<"CUDA error: "<<msg<<" : "<<cudaGetErrorString(e)<<"\n";
        std::exit(1);
    }
}
 
// Left and right ghosts
__global__ void bc_lr_kernel(GPUField U, int cur){          
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    // Ignore threads outside the valid row range
    if (j >= U.height) return;                              

    // Indices of the interior region
    int ng=U.ng;                                            
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;
    int jj=j;                                              

    // If j is within the range of the "ghost" rows, then use the nearest internal row
    if (jj<jb) jj=jb;                                       
    if (jj>jt) jj=jt;

    // Left ghost filling
    for(int i=0;i<ng;i++){                                  
        int id  = i + U.pitch*j;
        int src = il + U.pitch*jj;
        U.rho[cur][id]=U.rho[cur][src];
        U.rhou[cur][id]=U.rhou[cur][src];
        U.rhov[cur][id]=U.rhov[cur][src];
        U.E[cur][id]=U.E[cur][src];
    }
    
    // Right ghost filling
    for(int i=ir+1;i<=ir+ng;i++){                           
        int id  = i + U.pitch*j;
        int src = ir + U.pitch*jj;
        U.rho[cur][id]=U.rho[cur][src];
        U.rhou[cur][id]=U.rhou[cur][src];
        U.rhov[cur][id]=U.rhov[cur][src];
        U.E[cur][id]=U.E[cur][src];
    }
}

// Top and bottom ghosts
__global__ void bc_tb_kernel(GPUField U, int cur){          
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= U.pitch) return;

    int ng=U.ng;
    int il=ng, ir=ng+U.nx-1;
    int jb=ng, jt=ng+U.ny-1;

    int ii=i;
    if (ii<il) ii=il;
    if (ii>ir) ii=ir;

    for(int j=0;j<ng;j++){
        int id  = i + U.pitch*j;
        int src = ii + U.pitch*jb;
        U.rho[cur][id]=U.rho[cur][src];
        U.rhou[cur][id]=U.rhou[cur][src];
        U.rhov[cur][id]=U.rhov[cur][src];
        U.E[cur][id]=U.E[cur][src];
    }
    for(int j=jt+1;j<=jt+ng;j++){
        int id  = i + U.pitch*j;
        int src = ii + U.pitch*jt;
        U.rho[cur][id]=U.rho[cur][src];
        U.rhou[cur][id]=U.rhou[cur][src];
        U.rhov[cur][id]=U.rhov[cur][src];
        U.E[cur][id]=U.E[cur][src];
    }
}

void gpu_bc(GPUField& dU){
    
    // Create CUDA events for timing
    cudaEvent_t s,e;
    cudaEventCreate(&s); cudaEventCreate(&e);               

    // Grid configuration
    dim3 block(256);                                        
    dim3 grid_lr((dU.height + block.x - 1)/block.x);        
    dim3 grid_tb((dU.pitch  + block.x - 1)/block.x);        

    // CUDA event timing
    cudaEventRecord(s);                                     
    bc_lr_kernel<<<grid_lr, block>>>(dU, dU.cur);           
    ck(cudaGetLastError(), "bc_lr_kernel");                 
    bc_tb_kernel<<<grid_tb, block>>>(dU, dU.cur);           
    ck(cudaGetLastError(), "bc_tb_kernel");                 
    cudaEventRecord(e);                                     
    ck(cudaEventSynchronize(e), "sync bc");                 

    float ms=0.f;
    cudaEventElapsedTime(&ms, s, e);                        
    g_bc_ms += ms;

    // Destroy CUDA events
    cudaEventDestroy(s); cudaEventDestroy(e);               
}

// Return time
double gpu_get_bc_ms(){ return g_bc_ms; }                   
