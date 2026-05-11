#pragma once
#include "field.h"
#include "params.h"
#include "mpi_decomp.h"

// Compute dt for serial and OpenMP runs
double compute_dt_cfl(const Field& U, const Params& P);
// Compute dt for MPI and record the location of the minimum value
double compute_dt_cfl_mpi_strict(const Field& U, const Params& P, const MpiDecomp& D, int& i_min, int& j_min);        
