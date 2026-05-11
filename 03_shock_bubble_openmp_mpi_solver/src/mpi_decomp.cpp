#include "mpi_decomp.h"
#include <vector>
#include <cassert>
#include <iostream> 

// Compute starting index and local length
static inline void decomp_range(int N, int size, int rank, int& start, int& count) {
    int base = N / size;           // Minimum number of elements per rank
    int rem  = N % size;           // Remaining elements
    count = base + (rank < rem ? 1 : 0);                    // Extra element goes to the first rem ranks
    start = rank * base + (rank < rem ? rank : rem);        // Global start index for this rank
}

// MPI one-dimensional decomposition in the y-direction
MpiDecomp make_decomp_1d_y(int ny_global, int ng, MPI_Comm comm) {
    MpiDecomp D;
    MPI_Comm_rank(comm, &D.rank);
    MPI_Comm_size(comm, &D.size);
    D.ny_global = ny_global;            

    // Distribute ny_global evenly among size ranks
    decomp_range(ny_global, D.size, D.rank, D.j0, D.local_ny);

    // Neighboring ranks in the y-direction
    D.down = (D.rank == 0) ? MPI_PROC_NULL : D.rank - 1;
    D.up   = (D.rank == D.size - 1) ? MPI_PROC_NULL : D.rank + 1;
    return D;
}

// pack ng rows, but only interior columns (nx), for all 4 arrays
static void pack_rows_interior(const Field& U, int j_begin, int ng, std::vector<double>& buf) {
    int nx = U.nx;
    int ilo = U.ng;
    buf.resize((size_t)ng * nx * 4);
    size_t p = 0;
    // Outer loop: pack ng rows
    for (int g = 0; g < ng; ++g) {
        int j = j_begin + g;
        // Inner loop: Pack nx interior columns in each row
        for (int ii = 0; ii < nx; ++ii) {
            int i = ilo + ii;              // interior columns only
            int id = U.idx(i, j);
            buf[p++] = U.rho[id];
            buf[p++] = U.rhou[id];
            buf[p++] = U.rhov[id];
            buf[p++] = U.E[id];
        }
    }
}

// unpack ng rows, but only interior columns (nx), for all 4 arrays
static void unpack_rows_interior(Field& U, int j_begin, int ng, const std::vector<double>& buf) {
    int nx = U.nx;
    int ilo = U.ng;
    size_t p = 0;
    // Unpack data back into the field
    for (int g = 0; g < ng; ++g) {
        int j = j_begin + g;
        for (int ii = 0; ii < nx; ++ii) {
            int i = ilo + ii;
            int id = U.idx(i, j);
            U.rho[id]  = buf[p++];
            U.rhou[id] = buf[p++];
            U.rhov[id] = buf[p++];
            U.E[id]    = buf[p++];
        }
    }
}

// Exchange ghost rows with neighboring ranks
void exchange_halos_y(Field& U, const MpiDecomp& D, int ng, MPI_Comm comm) {

/*
    static int call = 0;
    if (call < 20) {        // Only print the first 20 times to avoid screen overflow
        int rank; MPI_Comm_rank(comm, &rank);
        if (rank == 0) std::cout << "[DBG] exchange_halos_y call=" << call << "\n";
    }
    call++;
*/

    // Internal area y index range (including ghost):
    // interior j: [ng, ng + local_ny - 1]
    // bottom ghost j: [0, ng-1]
    // top ghost    j: [ng+local_ny, ng+local_ny+ng-1]
    int j_interior_lo = U.ng;
    int j_interior_hi = U.ng + U.ny - 1;

    // Send bottom interior rows to the lower neighbor's top ghost region
    // Send top interior rows to the upper neighbor's bottom ghost region
    std::vector<double> send_down, send_up, recv_down, recv_up;

    // pack bottom interior rows
    pack_rows_interior(U, j_interior_lo, ng, send_down);
    // pack top interior rows
    pack_rows_interior(U, j_interior_hi - ng + 1, ng, send_up);

    recv_down.resize(send_down.size());
    recv_up.resize(send_up.size());

    MPI_Request reqs[4];
    int r = 0;

    // recv from down into bottom ghost
    MPI_Irecv(recv_down.data(), (int)recv_down.size(), MPI_DOUBLE, D.down, 100, comm, &reqs[r++]);
    // recv from up into top ghost
    MPI_Irecv(recv_up.data(), (int)recv_up.size(), MPI_DOUBLE, D.up,   101, comm, &reqs[r++]);

    // send to down (they receive as up)
    MPI_Isend(send_down.data(), (int)send_down.size(), MPI_DOUBLE, D.down, 101, comm, &reqs[r++]);
    // send to up (they receive as down)
    MPI_Isend(send_up.data(),   (int)send_up.size(),   MPI_DOUBLE, D.up,   100, comm, &reqs[r++]);

    MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);

    // unpack into ghost regions (only if neighbor exists)
    if (D.down != MPI_PROC_NULL) unpack_rows_interior(U, 0, ng, recv_down);
    if (D.up   != MPI_PROC_NULL) unpack_rows_interior(U, U.ng + U.ny, ng, recv_up);
    
}
