#pragma once
#include "field.h"
#include "params.h"
#include <string>

// CSV output for serial and OpenMP
// Write the density field in CSV format
void write_rho_csv(const Field& U, const std::string& path);                            
// Compute pressure from the conserved variables and write CSV
void write_p_csv  (const Field& U, const Params& P, const std::string& path);           

