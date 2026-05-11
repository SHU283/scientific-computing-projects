#include "field.h"
#include "params.h"
#include "init.h"
#include "bc.h"
#include "dt.h"
#include "hrsc.h"
#include "checks.h"
#include "timer.h"
#include "mpi_decomp.h"
#include "mpi_output.h"

#include <iostream>
#include <filesystem>
#include <mpi.h>
#include <algorithm>
#include <string>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    Params P;

    // weak scaling
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--weak" && i + 1 < argc) {
            int N = std::max(1, std::atoi(argv[i+1]));
            P.weakN = N;
            P.ny = 197 * N;
            P.Ly = 0.089 * N;
            ++i;
        }
    }

    // One rank per line
    MpiDecomp D = make_decomp_1d_y(P.ny, P.ng, comm);

    // rank0 prints some decomposition information
    if (D.rank == 0) {
        std::cout << "[MPI] size = " << D.size << "\n";
        std::cout << "[MPI] local_ny example (rank0) = " << D.local_ny << "\n";
    }

    // local field
    Field U(P.nx, D.local_ny, P.ng);

    // Create output directory (rank 0 only)
    if (D.rank == 0) std::filesystem::create_directories("output/scaling");

    // Initialize the local field from the global geometry
    init_shock_bubble_local(U, P, D.j0); 

    // Initial ghost synchronization
    apply_bc_lr_transmissive(U);
    exchange_halos_y(U, D, P.ng, comm);

    // Initialize time-stepping variables and timers
    double t = 0.0;
    int step = 0;
    int max_steps = 1000000;

    Timer t_dt, t_adv;

    // main time loop
    while (t < P.T_final && step < max_steps) {

        // dt: local then global min
        t_dt.start();
        double dt_local = compute_dt_cfl(U, P);
        double dt = 0.0;
        MPI_Allreduce(&dt_local, &dt, 1, MPI_DOUBLE, MPI_MIN, comm);
        t_dt.stop();

        if (t + dt > P.T_final) dt = P.T_final - t;

        // advance
        t_adv.start();
        advance_muscl_hancock_mpi(U, dt, P, D, comm);
        t_adv.stop();

        t += dt;
        step++;

        // Perform a numerical check every 200 steps
        if (step % 200 == 0) {
            int ok_local = check_physical(U, P, step, t) ? 1 : 0;
            int ok_global = 1;
            MPI_Allreduce(&ok_local, &ok_global, 1, MPI_INT, MPI_MIN, comm);

            if (D.rank == 0) {
                std::cout << "[MPI CHECK] step=" << step << " t=" << t
                          << " dt=" << dt << " ok=" << ok_global << "\n";
            }
            if (!ok_global) break;
        }
    }

    // output (rank0 writes)
    mpi_gather_and_write_csv(U, P, D, comm, "output/scaling/mpi_np" + std::to_string(D.size));
    //mpi_gather_and_write_csv(U, P, D, comm, "output/mpi_T2e-4_np" + std::to_string(D.size));

    // Rank0 Prints summary timing information
    if (D.rank == 0) {
        std::cout << "Done. steps=" << step << " final t=" << t << "\n";
        std::cout << "[TIMER MPI] dt   total = " << t_dt.acc  << " s\n";
        std::cout << "[TIMER MPI] adv  total = " << t_adv.acc << " s\n";
        std::cout << "[TIMER MPI] steps = " << step << "\n";
        std::cout << "[TIMER MPI] per-step adv = "
                  << (t_adv.acc / std::max(1, step)) << " s\n";
    }

    MPI_Finalize();
    return 0;
}
