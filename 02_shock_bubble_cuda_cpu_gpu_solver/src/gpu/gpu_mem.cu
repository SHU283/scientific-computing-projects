#include "core/gpu_field.h"
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <cstdlib>

// Check
static void ck(cudaError_t e, const char* msg) {
    if (e != cudaSuccess) {
        std::cerr << "CUDA error: " << msg << " : "
                  << cudaGetErrorString(e) << "\n";
        std::exit(1);
    }
}

// Allocate GPU memory for GPUField
void gpu_alloc(GPUField& dU, const Field& hU) {                                                  
    dU.nx = hU.nx; dU.ny = hU.ny; dU.ng = hU.ng;                                        
    dU.pitch = hU.pitch; dU.height = hU.height;
    dU.cur = 0;

    // Allocate two sets of arrays (double buffering) for each conserved quantity on the GPU
    size_t bytes = sizeof(double) * (size_t)gpu_size_tot(dU);                             
    for (int k=0;k<2;k++){                                                                 
        ck(cudaMalloc(&dU.rho[k],  bytes), "malloc rho");
        ck(cudaMalloc(&dU.rhou[k], bytes), "malloc rhou");
        ck(cudaMalloc(&dU.rhov[k], bytes), "malloc rhov");
        ck(cudaMalloc(&dU.E[k],    bytes), "malloc E");
    }
}

// Release the GPU memory and clear the structure
void gpu_free(GPUField& dU) {                                                                       
    for (int k=0;k<2;k++){
        cudaFree(dU.rho[k]);  cudaFree(dU.rhou[k]);
        cudaFree(dU.rhov[k]); cudaFree(dU.E[k]);
    }
    dU = {};
}

// CPU → GPU（Host to Device）
void gpu_copy_c2g(GPUField& dU, const Field& hU) {                                                
    int c = dU.cur;                                                                            
    size_t bytes = sizeof(double) * (size_t)gpu_size_tot(dU);                                     
    ck(cudaMemcpy(dU.rho[c],  hU.rho.data(),  bytes, cudaMemcpyHostToDevice), "h2d rho");        
    ck(cudaMemcpy(dU.rhou[c], hU.rhou.data(), bytes, cudaMemcpyHostToDevice), "h2d rhou");
    ck(cudaMemcpy(dU.rhov[c], hU.rhov.data(), bytes, cudaMemcpyHostToDevice), "h2d rhov");
    ck(cudaMemcpy(dU.E[c],    hU.E.data(),    bytes, cudaMemcpyHostToDevice), "h2d E");
}

// GPU → CPU（Device to Host）
void gpu_copy_g2c(Field& hU, const GPUField& dU) {                                             
    int c = dU.cur;
    size_t bytes = sizeof(double) * (size_t)gpu_size_tot(dU);
    ck(cudaMemcpy(hU.rho.data(),  dU.rho[c],  bytes, cudaMemcpyDeviceToHost), "d2h rho");
    ck(cudaMemcpy(hU.rhou.data(), dU.rhou[c], bytes, cudaMemcpyDeviceToHost), "d2h rhou");
    ck(cudaMemcpy(hU.rhov.data(), dU.rhov[c], bytes, cudaMemcpyDeviceToHost), "d2h rhov");
    ck(cudaMemcpy(hU.E.data(),    dU.E[c],    bytes, cudaMemcpyDeviceToHost), "d2h E");
}
