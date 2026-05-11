#include "mpi_output.h"
#include "eos.h"          // pressure_from_cons
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

// Pack the interior region (excluding ghosts) into contiguous buffers of size nx*ny_local
static void pack_interior(const Field& U, std::vector<double>& rho,
                          std::vector<double>& rhou,
                          std::vector<double>& rhov,
                          std::vector<double>& E)
{
    int nx = U.nx;
    int ny = U.ny;
    int ng = U.ng;

    // Allocate output vectors
    rho.resize((size_t)nx * ny);
    rhou.resize((size_t)nx * ny);
    rhov.resize((size_t)nx * ny);
    E.resize((size_t)nx * ny);

    // Flatten the 2D interior region into 1D buffers
    size_t p = 0;
    for (int j = 0; j < ny; ++j) {
        int jj = ng + j;
        for (int i = 0; i < nx; ++i) {
            int ii = ng + i;
            int id = U.idx(ii, jj);
            rho[p]  = U.rho[id];
            rhou[p] = U.rhou[id];
            rhov[p] = U.rhov[id];
            E[p]    = U.E[id];
            ++p;
        }
    }
}

// Write a 1D array to CSV using a 2D grid layout
static void write_csv_2d(const std::vector<double>& a, int nx, int ny, const std::string& path)
{
    std::ofstream f(path);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            f << a[(size_t)j * nx + i];
            if (i != nx - 1) f << ",";
        }
        f << "\n";
    }
}

// Collect the local subdomains (excluding the ghost regions) of each MPI rank (in the interior)
// gather them to rank 0 to assemble the global field
// Calculate the pressure on rank 0 and write the global rho/p in CSV format
void mpi_gather_and_write_csv(const Field& U_local,
                              const Params& P,
                              const MpiDecomp& D,
                              MPI_Comm comm,
                              const std::string& out_prefix)
{
    // Pack local interior data
    std::vector<double> lrho, lrhou, lrhov, lE;
    pack_interior(U_local, lrho, lrhou, lrhov, lE);

    int nx = P.nx;
    int ny = P.ny;

    // counts/displs for Gatherv (units: doubles)
    std::vector<int> counts, displs;
    if (D.rank == 0) {
        counts.resize(D.size);
        displs.resize(D.size);
        for (int r = 0; r < D.size; ++r) {
            // same decomposition rule as make_decomp_1d_y()
            int base = ny / D.size;
            int rem  = ny % D.size;
            int local_ny = base + (r < rem ? 1 : 0);
            int j0 = r * base + (r < rem ? r : rem);

            counts[r] = nx * local_ny;
            displs[r] = nx * j0;
        }
    }

    // root buffers (global conserved)
    std::vector<double> grho, grhou, grhov, gE;
    if (D.rank == 0) {
        grho.resize((size_t)nx * ny);
        grhou.resize((size_t)nx * ny);
        grhov.resize((size_t)nx * ny);
        gE.resize((size_t)nx * ny);
    }

    // MPI_Gatherv for each of the four variables
    MPI_Gatherv(lrho.data(),  (int)lrho.size(),  MPI_DOUBLE,
                D.rank == 0 ? grho.data()  : nullptr,
                D.rank == 0 ? counts.data(): nullptr,
                D.rank == 0 ? displs.data(): nullptr,
                MPI_DOUBLE, 0, comm);

    MPI_Gatherv(lrhou.data(), (int)lrhou.size(), MPI_DOUBLE,
                D.rank == 0 ? grhou.data() : nullptr,
                D.rank == 0 ? counts.data(): nullptr,
                D.rank == 0 ? displs.data(): nullptr,
                MPI_DOUBLE, 0, comm);

    MPI_Gatherv(lrhov.data(), (int)lrhov.size(), MPI_DOUBLE,
                D.rank == 0 ? grhov.data() : nullptr,
                D.rank == 0 ? counts.data(): nullptr,
                D.rank == 0 ? displs.data(): nullptr,
                MPI_DOUBLE, 0, comm);

    MPI_Gatherv(lE.data(),    (int)lE.size(),    MPI_DOUBLE,
                D.rank == 0 ? gE.data()    : nullptr,
                D.rank == 0 ? counts.data(): nullptr,
                D.rank == 0 ? displs.data(): nullptr,
                MPI_DOUBLE, 0, comm);

    if (D.rank != 0) return;

    // Reconstruct pressure on rank 0
    std::vector<double> gp((size_t)nx * ny);
    for (size_t k = 0; k < gp.size(); ++k) {
        double rho  = grho[k];
        double rhou = grhou[k];
        double rhov = grhov[k];
        double E    = gE[k];
        gp[k] = pressure_from_cons(rho, rhou, rhov, E, P.gamma);
    }

    // Write global CSV files
    write_csv_2d(grho, nx, ny, out_prefix + "_rho.csv");
    write_csv_2d(gp,   nx, ny, out_prefix + "_p.csv");
}
