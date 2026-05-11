#include "hrsc.h"
#include "bc.h"
#include "mpi_decomp.h"

// Copy the values of row src_j into ghost rows in the y-direction
// dir = -1: down to (jb - g)
// dir = +1: up to (jt + g)
static inline void copy_row_to_ghosts(Field& U, int src_j, int dir)
{
    int ng = U.ng;
    for (int i = 0; i < U.pitch; ++i) {

        // The value of src remains constant within the g loop
        int src = U.idx(i, src_j);

        for (int g = 1; g <= ng; ++g) {
            int dst = U.idx(i, src_j + dir * g);

            U.rho[dst]  = U.rho[src];
            U.rhou[dst] = U.rhou[src];
            U.rhov[dst] = U.rhov[src];
            U.E[dst]    = U.E[src];
        }
    }
}

// Fill top and bottom ghost rows at the global physical boundaries
static inline void apply_bc_tb_global(Field& U, const MpiDecomp& D)
{
    int ng = U.ng;
    int jb = ng;             // bottom internal boundary j
    int jt = ng + U.ny - 1;  // top internal boundary j

    // Bottom physical boundary: only handled by the bottom rank
    if (D.down == MPI_PROC_NULL) {
        copy_row_to_ghosts(U, jb, -1);
    }

    // Top physical boundary: only handled by the top rank
    if (D.up == MPI_PROC_NULL) {
        copy_row_to_ghosts(U, jt, +1);
    }
}

// Update ghost cells from neighbors and physical boundaries
static inline void sync_ghosts(Field& U, const MpiDecomp& D, MPI_Comm comm)
{
    apply_bc_lr_transmissive(U);        // Fill global left/right first.
    exchange_halos_y(U, D, U.ng, comm); // Exchange top/bottom ghost rows with neighboring ranks
    apply_bc_tb_global(U, D);           // Fill global top/bottom(Only the boundary rank)
}

// MPI implementation of the second-order method
void advance_muscl_hancock_mpi(Field& U, double dt, const Params& P,
                              const MpiDecomp& D, MPI_Comm comm)
{
    double dx = P.Lx / P.nx;
    double dy = P.Ly / P.ny;

    double half_dt = 0.5 * dt;

    // Strang splitting, x(half) -> y(full) -> x(half)
    sync_ghosts(U, D, comm);
    sweep_x(U, half_dt, dx, P.gamma);

    sync_ghosts(U, D, comm);
    sweep_y(U, dt, dy, P.gamma);

    sync_ghosts(U, D, comm);
    sweep_x(U, half_dt, dx, P.gamma);

    sync_ghosts(U, D, comm);
}
