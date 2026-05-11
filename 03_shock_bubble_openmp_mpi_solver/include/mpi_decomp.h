#pragma once
#include "field.h"
#include <mpi.h>

// MPI domain decomposition
struct MpiDecomp {
    int rank = 0, size = 1;
    int ny_global = 0;          // Global number of grid points in the y direction
    int j0 = 0;                 // Global starting y-index of this rank, excluding ghosts
    int local_ny = 0;           // Local y-size of this rank, excluding ghosts
    int up = MPI_PROC_NULL;     // Upper neighboring rank
    int down = MPI_PROC_NULL;   // Lower neighboring rank
};

// Build the 1D y-direction decomposition for each rank
MpiDecomp make_decomp_1d_y(int ny_global, int ng, MPI_Comm comm);               

// Exchange ghost rows with the upper and lower ranks
void exchange_halos_y(Field& U, const MpiDecomp& D, int ng, MPI_Comm comm);     