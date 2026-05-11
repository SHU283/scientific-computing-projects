#pragma once
#include "field.h"

// Split BC for MPI usage
void apply_bc_lr_transmissive(Field& U);        // Set the left and right ghost points
void apply_bc_tb_transmissive(Field& U);        // Set the top and bottom ghost points

// serial/omp BC
void apply_bc_transmissive(Field& U);              
