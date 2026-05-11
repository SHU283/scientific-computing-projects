#pragma once
#include "field.h"
#include "params.h"

// Initialize the calculation domain (serial / OpenMP)
void init_shock_bubble(Field& U, const Params& P);                              

// Initialize the local domain for MPI
void init_shock_bubble_local(Field& U, const Params& P, int j0_global);         

