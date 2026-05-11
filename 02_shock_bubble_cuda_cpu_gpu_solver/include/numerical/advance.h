#pragma once
#include "core/field.h"

// First-order Godunov (HLL flux) + Explicit Euler time advancement
void advance_1st(Field& U, double dt, double dx, double dy, double gamma);

// second order：MUSCL-Hancock + Strang splitting (x half, y full, x half)
void advance_2nd(Field& U, double dt, double dx, double dy, double gamma);