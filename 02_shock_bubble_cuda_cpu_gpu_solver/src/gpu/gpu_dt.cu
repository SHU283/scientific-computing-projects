#include "core/gpu_field.h"
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

// Check
static void ck(cudaError_t e, const char* msg){                                         
    if(e!=cudaSuccess){
        std::cerr<<"CUDA error: "<<msg<<" : "<<cudaGetErrorString(e)<<"\n";
        std::exit(1);
    }
}

// Calculate p using E
__device__ __forceinline__                                                              
double E_to_p_d(double rho, double rhou, double rhov, double E, double g){
    double u=rhou/rho, v=rhov/rho;
    return (g-1.0)*(E - 0.5*rho*(u*u+v*v));
}

// Calculate the max safe dt for a single cell
__device__ __forceinline__                                                              
double safe_dt_cell(double rho,double rhou,double rhov,double E, double dx,double dy,double g){
    // Check rho and p
    if(rho<=0.0) return 1e30;                                                         
    double u=rhou/rho, v=rhov/rho;
    double p=E_to_p_d(rho,rhou,rhov,E,g);
    if(p<=0.0) return 1e30;
    // Sound speed c
    double c = sqrt(fmax(0.0, g*p/rho));
    // Maximum wave speed in x and y
    double smax_x = fabs(u)+c;                                                         
    double smax_y = fabs(v)+c;
    if(smax_x<=0.0 || smax_y<=0.0) return 1e30;
    // Return the upper limit of the dt of this cell
    // Still need to multiply CFL number
    return fmin(dx/smax_x, dy/smax_y);                                                  
}

// A round of minimum value reduction
// For each block, take the minimum value
// Finally, each block outputs the minimum dt to partial[blockIdx.x]
__global__ void dt_partial_min_kernel(GPUField U, int cur, double dx,               
                                        double dy,double g, double* partial){           // partial stores dt for each block     
    // Dynamic shared memory array, shared data within blocks
    extern __shared__ double s[];     
    // Thread-local and global indices                                               
    int tid = threadIdx.x;                                                              
    int idx = blockIdx.x * blockDim.x + tid;                               

    // Interior only
    int N = U.nx * U.ny;      
    // Safe initial dt                                                        
    double src = 1e30;                                                                 

    // Map the one-dimensional idx back to the interior's (ii, jj), and then map it to the ghosted (i, j)
    if(idx < N){                                                                        
        int ii = idx % U.nx;                                                            
        int jj = idx / U.nx;                                                           
        int i = U.ng + ii;                                                           
        int j = U.ng + jj;
        int id = i + U.pitch*j;    
        // Calculate dt                                                
        src = safe_dt_cell(U.rho[cur][id],U.rhou[cur][id],U.rhov[cur][id],U.E[cur][id], dx,dy,g);    
    }

    // Store each thread's dt in shared memory
    s[tid]=src;                                                                    
    __syncthreads();   

    // block internal reduction (halving and taking the minimum)
    for(int s2=blockDim.x/2; s2>0; s2>>=1){                                             
        if(tid < s2) s[tid] = fmin(s[tid], s[tid+s2]);
        __syncthreads();
    }
    // Write the minimum value for each block to partial
    if(tid==0) partial[blockIdx.x]=s[0];                                                
}

// Global minimum value reduction kernel
__global__ void reduce_min_kernel(const double* in,double* out,int n){             
    // Shared memory: Temporary storage within a block     
    extern __shared__ double s[];
    // Each thread takes one input element (if it exceeds the limit, it assigns 1e30)
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + tid;
    double v = (idx < n) ? in[idx] : 1e30;                                             
    // Put the value of each thread into the shared variable and synchronize
    s[tid]=v;                                                                      
    __syncthreads();
    // The folding method is used to find the smallest dt in the current block
    for(int s2=blockDim.x/2; s2>0; s2>>=1){                                             
        if(tid < s2) s[tid] = fmin(s[tid], s[tid+s2]); 
        __syncthreads();
    }
    // Assign the smallest value of each block to the output array
    if(tid==0) out[blockIdx.x]=s[0];                                                    
}

// Calculate the global minimum stable time step dt_min on the GPU, and then multiply it by CFL to obtain the result
double gpu_compute_dt_cfl(const GPUField& dU,double dx,double dy,double gamma,double CFL){         
    // Calculate the total number of cells + configure blocks/grid
    int N = dU.nx * dU.ny;                                                                          
    int block = 256;                                                                                
    int grid  = (N + block - 1) / block;                                                            

    // Allocate two temporary arrays, bufA and bufB, on the GPU to separate the input and output
    double *bufA=nullptr, *bufB=nullptr;                                                            
    ck(cudaMalloc(&bufA, sizeof(double)*grid), "malloc bufA");                                    
    ck(cudaMalloc(&bufB, sizeof(double)*grid), "malloc bufB");        

    // Calculate dt for the entire grid and output the minimum value of each block to bufA
    dt_partial_min_kernel<<<grid,block, sizeof(double)*block>>>(dU, dU.cur, dx, dy, gamma, bufA);   
    ck(cudaGetLastError(), "dt_partial_min_kernel");                                                

    // Multiple rounds of reduction: Reduce bufA to contain only one value
    int curN = grid;                                                                              
    bool flip = false;                                                                             
    while(curN > 1){                                                                           
        int g2 = (curN + block - 1) / block;                                                   
        const double* in  = flip ? bufB : bufA;                                         
        double*       out = flip ? bufA : bufB;
        reduce_min_kernel<<<g2,block,sizeof(double)*block>>>(in, out, curN);                  
        ck(cudaGetLastError(), "reduce_min_kernel");
        curN = g2;                                                                             
        flip = !flip;                                                                    
    }

    // Copy the last dt_min from the GPU back to the CPU
    double dt_min=1e30;                                                                      
    const double* last = flip ? bufB : bufA;
    ck(cudaMemcpy(&dt_min, last, sizeof(double), cudaMemcpyDeviceToHost), "copy dt_min");

    // Release temporary GPU memory
    cudaFree(bufA); cudaFree(bufB);                                                          

    // Fallback if no valid minimum time step was found
    if(!(dt_min < 1e29)) dt_min = 1e-12;                                                  
    return CFL * dt_min;
}
