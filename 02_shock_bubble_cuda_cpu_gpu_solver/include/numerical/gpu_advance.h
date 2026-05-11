#pragma once
#include "core/gpu_field.h"

// 1st order (Forward Euler)
void gpu_advance_1st(GPUField& dU, double dt, double dx, double dy, double gamma);  

// 2nd order (MUSCL-Hancock)
void gpu_advance_2nd(GPUField& dU, double dt, double dx, double dy, double gamma);  
