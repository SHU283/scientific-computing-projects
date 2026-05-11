#include "core/field.h"
#include <iostream>
#include <filesystem> 
#include "numerical/advance.h" 
#include <chrono>
#include "numerical/bc.h"
#include "numerical/dt.h"
#include "physics/init.h"
#include "io/output.h"


int main(){
    // Grid size and field variables
    int nx = 500, ny = 197;     
    Field U(nx, ny, 2);       

    // Physical parameter
    double gamma = 1.4;       
    double p0 = 1.01325e5;    
    double rho_air = 1.29, rho_he = 0.214;     

    // Geometrical parameter
    double Lx = 0.225, Ly = 0.089;
    double xb = 0.035, yb = 0.0, rb = 0.025;
    double xshock = 0.005;

    // Initialize the field
    init(U, Lx, Ly, gamma, p0, rho_air, rho_he, xb, yb, rb, xshock);

    // Grid step size and CFL
    double dx = Lx / nx;
    double dy = Ly / ny;
    double CFL = 0.4;

    // Variables related to time progression
    double t = 0.0;
    double tEnd = 1.0e-4;
    int outEvery = 10;   // set a large number for benchmark
    int step = 0;

    bool use_second_order = false;   // Second order = true, first order = false
    // Write generated fields under a local data directory so the project can
    // be moved or submitted without editing machine-specific paths.
    const std::filesystem::path outBase = "data";
    const std::filesystem::path outDir = outBase / (use_second_order ? "data_2nd" : "data_1st");
    std::filesystem::create_directories(outDir);

    
    // Output the initial frame
    write_rho_csv(U, (outDir / "rho_0.csv").string());
    write_p_csv(U, (outDir / "p_0.csv").string(), gamma);

    // Timing variable
    auto t0 = std::chrono::high_resolution_clock::now();  
    double t_bc = 0.0, t_dt = 0.0, t_step = 0.0;  

    while (t < tEnd) {

        // Run BC and calculate the time
        {
            auto a = std::chrono::high_resolution_clock::now();
            fill_ghost(U);
            auto b = std::chrono::high_resolution_clock::now();
            t_bc += std::chrono::duration<double>(b - a).count();
        }

        // Run dt and calculate the time
        double dt;
        {
            auto a = std::chrono::high_resolution_clock::now();
            dt = cfl_dt(U, dx, dy, gamma, CFL);
            auto b = std::chrono::high_resolution_clock::now();
            t_dt += std::chrono::duration<double>(b - a).count();
        }

        if (t + dt > tEnd) dt = tEnd - t;

        // Run one step and calculate the time
        {
            auto a = std::chrono::high_resolution_clock::now();
            if (use_second_order) {
                advance_2nd(U, dt, dx, dy, gamma);
            } else {
                advance_1st(U, dt, dx, dy, gamma);
            }
            auto b = std::chrono::high_resolution_clock::now();
            t_step += std::chrono::duration<double>(b - a).count();
        }

        t += dt;
        step++;

        // Output every outEvery steps
        if (step % outEvery == 0) {
            write_rho_csv(U, (outDir / ("rho_" + std::to_string(step) + ".csv")).string());
            write_p_csv(U, (outDir / ("p_" + std::to_string(step) + ".csv")).string(), gamma);
            std::cout << "step=" << step << " t=" << t << " dt=" << dt << "\n";
        }
    }

    // Calculate the total time and print the result
    auto t1 = std::chrono::high_resolution_clock::now();
    double seconds = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "[TIMER] CPU loop time = " << seconds << " s\n";
    std::cout << "[TIMER] CPU bc total   = " << t_bc   << " s\n";
    std::cout << "[TIMER] CPU dt total   = " << t_dt   << " s\n";
    std::cout << "[TIMER] CPU step total = " << t_step << " s\n";

    // Output every outEvery steps
    std::cout << "Done. steps=" << step << " final t=" << t << "\n";
    return 0;
}
 
