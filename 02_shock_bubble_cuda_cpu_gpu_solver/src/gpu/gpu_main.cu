#include <iostream>
#include <filesystem>
#include <chrono>
#include <cuda_runtime.h>          

#include "core/field.h"
#include "core/gpu_field.h"
#include "numerical/gpu_advance.h"
#include "physics/init.h"          
#include "io/output.h"             


int main(){
    // Set up the CPU grid and initial field
    int nx = 500, ny = 197;                                                              
    Field U(nx, ny, 2);

    // Physical constants + benchmark parameters
    double gamma = 1.4;  
    double p0 = 1.01325e5;
    double rho_air = 1.29, rho_he = 0.214;

    double Lx = 0.225, Ly = 0.089;
    double xb = 0.035, yb = 0.0, rb = 0.025;
    double xshock = 0.005;

    // Initialization (cpu)
    init(U, Lx, Ly, gamma, p0, rho_air, rho_he, xb, yb, rb, xshock);                      

    // Grid step size and CFL parameter
    double dx = Lx / nx;                                                            
    double dy = Ly / ny;
    double CFL = 0.4;

    // Time cycle control parameters
    double t = 0.0;                                                                      
    double tEnd = 1.0e-4;
    int outEvery = 10; 
    int step = 0;

    // Select first/second order + output directory
    bool use_second_order = false;                                                        
    // Write generated fields under a local data directory so the project can
    // be moved or submitted without editing machine-specific paths.
    const std::filesystem::path outBase = "data";
    const std::filesystem::path outDir = outBase / (use_second_order ? "data_gpu_2nd" : "data_gpu_1st");
    std::filesystem::create_directories(outDir);

    // Allocate the GPU field and copy the initial CPU data to the GPU
    GPUField dU{};                                                                      
    gpu_alloc(dU, U);
    gpu_copy_c2g(dU, U);

    // Three types of timers: dt / bc / step (CPU wall time)
    double t_dt   = 0.0;                                                                    
    double t_bc_wall   = 0.0;                                                              
    double t_step_wall = 0.0;                                                             

    // initial output (host)                                                                
    write_rho_csv(U, (outDir / "rho_0.csv").string());
    write_p_csv(U, (outDir / "p_0.csv").string(), gamma);

    // Start total wall-time measurement for the main loop
    auto t0 = std::chrono::high_resolution_clock::now();                               

    // Main loop: BC → dt → step → output
    while (t < tEnd) {                                                                  

        // BC
        {
            auto a = std::chrono::high_resolution_clock::now();                           
            gpu_bc(dU);
            cudaDeviceSynchronize();                                                    
            auto b = std::chrono::high_resolution_clock::now();
            t_bc_wall += std::chrono::duration<double>(b - a).count();
        }

        // dt
        double dt;
        {
            auto a = std::chrono::high_resolution_clock::now();                              
            dt = gpu_compute_dt_cfl(dU, dx, dy, gamma, CFL);
            cudaDeviceSynchronize(); 
            auto b = std::chrono::high_resolution_clock::now();
            t_dt += std::chrono::duration<double>(b - a).count();
        }

        // Make sure not to exceed the final time
        if (t + dt > tEnd) dt = tEnd - t;                                                      

        // Advance step (1st or 2nd)
        {
            auto a = std::chrono::high_resolution_clock::now();                              
            if (use_second_order) {
                gpu_advance_2nd(dU, dt, dx, dy, gamma);
            } else {
                gpu_advance_1st(dU, dt, dx, dy, gamma);
            }
            cudaDeviceSynchronize();
            auto b = std::chrono::high_resolution_clock::now();
            t_step_wall += std::chrono::duration<double>(b - a).count();
        }

        // Update time and steps
        t += dt;                                                                                                   
        step++;

        // Write output every outEvery steps
        if (step % outEvery == 0) {
            gpu_copy_g2c(U, dU);                                                           
            write_rho_csv(U, (outDir / ("rho_" + std::to_string(step) + ".csv")).string());
            write_p_csv(U, (outDir / ("p_" + std::to_string(step) + ".csv")).string(), gamma);
            std::cout << "GPU step=" << step << " t=" << t << " dt=" << dt << "\n";
        }
    }

    // Pure GPU timing
    std::cout << "[TIMER] GPU bc total   = " << gpu_get_bc_ms()   * 1e-3 << " s\n";            
    std::cout << "[TIMER] GPU step total = " << gpu_get_step_ms() * 1e-3 << " s\n";         
    // CPU wall time
    std::cout << "[TIMER] GPU dt total   = " << t_dt << " s\n";                          
    // Step count statistics
    std::cout << "[TIMER] GPU step calls = " << gpu_get_step_calls() << "\n";            
    
    // Total wall time + Sub-item wall time
    auto t1 = std::chrono::high_resolution_clock::now();
    double seconds = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[TIMER] GPU total wall time = " << seconds << " s\n";                       
    std::cout << "[TIMER] GPU bc wall (sync)   = " << t_bc_wall   << " s\n";                
    std::cout << "[TIMER] GPU step wall (sync) = " << t_step_wall << " s\n";                    

    // Clear resources and terminate
    gpu_free(dU);                                                                   
    std::cout << "Done. steps=" << step << " final t=" << t << "\n";
    return 0;
}  
