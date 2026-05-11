#include "field.h"
#include "params.h"
#include "init.h"
#include "bc.h"
#include "dt.h"
#include "hrsc.h"
#include "output.h"
#include <iostream>
#include <filesystem>
#include "checks.h"
#include "timer.h"


int main()
{
    // Create parameters, grid, and output folder
    Params P;
    Field U(P.nx, P.ny, P.ng);

    std::filesystem::create_directories("output");

    // Initialize the field, fill ghost cells once, and write the initial output
    init_shock_bubble(U, P);
    apply_bc_transmissive(U);

    write_rho_csv(U, "output/rho_0.csv");
    write_p_csv(U, P, "output/p_0.csv");

    // Time loop timer
    double t = 0.0;
    int step = 0;
    //int max_steps = 2000;
    int max_steps = 1000000; 

    Timer t_bc, t_dt, t_adv;
    int steps_compute = 0;

    // Main loop: BC -> dt -> advance
    while (t < P.T_final && step < max_steps) {

        // BC
        t_bc.start();
        apply_bc_transmissive(U);
        t_bc.stop();

        // dt
        t_dt.start();
        double dt = compute_dt_cfl(U, P);
        t_dt.stop();

        if (t + dt > P.T_final) dt = P.T_final - t;

        // ADVANCE
        t_adv.start();
        advance_muscl_hancock(U, dt, P);
        t_adv.stop();

        t += dt;
        step++;
        steps_compute++;

        // Progress Printing + Stability check
        if (step % P.output_every == 0) {
            //write_rho_csv(U, "output/rho_" + std::to_string(step) + ".csv");
            //write_p_csv  (U, P, "output/p_"   + std::to_string(step) + ".csv");
            std::cout << "step=" << step << " t=" << t << " dt=" << dt << "\n";
        }

        if (step % 200 == 0) {
            if (!check_physical(U, P, step, t)) {
                std::cerr << "Physical check failed. Stop.\n";
                break;
            }
        }

    }

    // Output the final field
    write_rho_csv(U, "output/serial_rho_final.csv");
    write_p_csv(U, P, "output/serial_p_final.csv");

    // Print timing statistics
    std::cout << "Done. steps=" << step << " final t=" << t << "\n";
    std::cout << "\n[TIMER] bc   total = " << t_bc.acc << " s\n";
    std::cout << "[TIMER] dt   total = " << t_dt.acc << " s\n";
    std::cout << "[TIMER] adv  total = " << t_adv.acc << " s\n";
    std::cout << "[TIMER] steps = " << steps_compute << "\n";
    std::cout << "[TIMER] per-step adv = "
              << (t_adv.acc / std::max(1, steps_compute)) << " s\n";
    
    return 0;
}
