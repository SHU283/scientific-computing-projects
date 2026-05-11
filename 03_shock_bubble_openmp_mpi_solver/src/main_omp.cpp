#include "field.h"
#include "params.h"
#include "init.h"
#include "bc.h"
#include "dt.h"
#include "hrsc.h"
#include "output.h"
#include "checks.h"
#include "timer.h"

#include <iostream>
#include <filesystem>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    Params P;

    // Weak scaling: Controlled by environment variables
    if (const char* s = std::getenv("HRSC_WEAK_N")) {
        int N = std::max(1, std::atoi(s));
        P.weakN = N;
        P.ny = 197 * N;
        P.Ly = 0.089 * N;
    }

    // Allocate the entire grid Field (global single block, not decomposed)
    Field U(P.nx, P.ny, P.ng);

    // Create output directory
    std::filesystem::create_directories("output/omp");

// Print OpenMP thread information
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
        std::cout << "[OMP] using threads = " << omp_get_num_threads() << "\n";
    }
#endif

#ifdef _OPENMP
    std::cout << "[OMP] max threads = " << omp_get_max_threads() << "\n";
#endif

    // Initialize the shock-bubble field
    init_shock_bubble(U, P);

    //write_rho_csv(U, "output/rho_0.csv");
    //write_p_csv(U, P, "output/p_0.csv");

    // Time loop until T_final is reached
    double t = 0.0;
    int step = 0;
    //int max_steps = 2000;
    int max_steps = 1000000;   // Safety cap, unlikely to be reached

    Timer t_bc, t_dt, t_adv;
    int steps_compute = 0;

    while (t < P.T_final && step < max_steps) {

        // Refresh ghost cells with boundary conditions
        t_bc.start();
        apply_bc_transmissive(U);
        t_bc.stop();

        // Calculate the time step dt of the CFL condition
        t_dt.start();
        double dt = compute_dt_cfl(U, P);
        t_dt.stop();

        // Prevent stepping past the final time
        if (t + dt > P.T_final) dt = P.T_final - t;

        // Advance one step (core calculation)
        t_adv.start();
        advance_muscl_hancock(U, dt, P);
        t_adv.stop();

        // Update time / Count
        t += dt;
        step++;
        steps_compute++;

        // Print every 50 steps and check every 200 steps
        if (step % 50 == 0) {
            std::cout << "step=" << step << " t=" << t << " dt=" << dt << "\n";
        }


        if (step % 200 == 0) {
            if (!check_physical(U, P, step, t)) {
                std::cerr << "Physical check failed. Stop.\n";
                break;
            }
        }
    }

    // Actual number of threads used
    int nth = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
    nth = omp_get_num_threads();
    }
#endif

    // Output the final field
    write_rho_csv(U, "output/omp/rho_final_ompN" + std::to_string(nth) + ".csv");
    write_p_csv  (U, P, "output/omp/p_final_ompN" + std::to_string(nth) + ".csv");

    // Printed timing statistics
    std::cout << "Done. steps=" << step << " final t=" << t << "\n";
    std::cout << "\n[TIMER] bc   total = " << t_bc.acc << " s\n";
    std::cout << "[TIMER] dt   total = " << t_dt.acc << " s\n";
    std::cout << "[TIMER] adv  total = " << t_adv.acc << " s\n";
    std::cout << "[TIMER] steps = " << steps_compute << "\n";
    std::cout << "[TIMER] per-step adv = "
              << (t_adv.acc / std::max(1, steps_compute)) << " s\n";

    return 0;
}
