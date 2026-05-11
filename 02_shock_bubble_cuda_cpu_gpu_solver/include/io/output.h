#pragma once
#include <string>
#include "core/field.h"

// Write the pressure data into a CSV file
void write_p_csv(const Field& U, const std::string& path, double gamma);
// Write the density data into a CSV file
void write_rho_csv(const Field&, const std::string&);  
