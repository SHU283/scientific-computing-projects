#include "field.h"

Field::Field(int nx_, int ny_, int ng_) {
    nx = nx_;
    ny = ny_;
    ng = ng_;

    // Storage size
    pitch  = nx + 2*ng;
    height = ny + 2*ng;

    // Allocate and initialize four arrays of conserved quantities
    int N = pitch * height;
    rho.assign(N, 0.0);
    rhou.assign(N, 0.0);
    rhov.assign(N, 0.0);
    E.assign(N, 0.0);
}
