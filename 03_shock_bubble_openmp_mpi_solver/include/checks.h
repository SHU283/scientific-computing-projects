#pragma once
#include "field.h"
#include "params.h"

// Check for physical legality
bool check_physical(const Field& U, const Params& P, int step, double t);       
