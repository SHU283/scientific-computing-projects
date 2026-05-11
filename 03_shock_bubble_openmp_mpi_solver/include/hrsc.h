#pragma once
#include "field.h"
#include "params.h"
#include "mpi_decomp.h"
#include <mpi.h>

// Full advancement routine for serial/OpenMP
void advance_muscl_hancock(Field& U, double dt, const Params& P);

// Split sweeps for MPI, with exchange and boundary conditions inserted between them
void sweep_x(Field& U, double dt, double dx, double gamma);
void sweep_y(Field& U, double dt, double dy, double gamma);
void advance_muscl_hancock_mpi(Field& U, double dt, const Params& P,
                               const MpiDecomp& D, MPI_Comm comm);