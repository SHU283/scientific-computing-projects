#pragma once
#include "field.h"
#include "params.h"
#include "mpi_decomp.h"
#include <mpi.h>
#include <string>

// MPI data writing to file
void mpi_gather_and_write_csv(const Field& U_local,                     
                              const Params& P,                          
                              const MpiDecomp& D,                       
                              MPI_Comm comm,
                              const std::string& out_prefix); 

