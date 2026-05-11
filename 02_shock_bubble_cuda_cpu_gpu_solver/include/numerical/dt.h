#pragma once
#include "core/field.h"

// Compute the CFL time step
double cfl_dt(const Field& U, double dx, double dy, double gamma, double CFL);   
