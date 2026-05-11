#include "core/field.h"
#include <algorithm>

void fill_ghost(Field& U)
{
    const int ng = U.ng;
    const int nx = U.nx;
    const int ny = U.ny;

    // Indices of the physical boundary cells
    const int left_edge   = ng;
    const int right_edge  = ng + nx - 1;
    const int bottom_edge = ng;
    const int top_edge    = ng + ny - 1;

    // Fill ghost columns on the left and right boundaries
    for (int j = 0; j < U.height; ++j) {
        for (int i = 0; i < ng; ++i) {
            // Left boundary: copy the nearest interior cell
            int ghost_L = U.index(i, j);
            int valid_L = U.index(left_edge, j);
            U.rho[ghost_L]  = U.rho[valid_L];
            U.rhou[ghost_L] = U.rhou[valid_L];
            U.rhov[ghost_L] = U.rhov[valid_L];
            U.E[ghost_L]    = U.E[valid_L];

            // Right boundary: copy the nearest interior cell
            int ghost_R = U.index(right_edge + 1 + i, j);
            int valid_R = U.index(right_edge, j);
            U.rho[ghost_R]  = U.rho[valid_R];
            U.rhou[ghost_R] = U.rhou[valid_R];
            U.rhov[ghost_R] = U.rhov[valid_R];
            U.E[ghost_R]    = U.E[valid_R];
        }
    }

    // Fill ghost rows on the bottom and top boundaries
    for (int i = 0; i < U.pitch; ++i) {
        for (int j = 0; j < ng; ++j) {
            // Bottom boundary: copy the nearest interior cell
            int ghost_B = U.index(i, j);
            int valid_B = U.index(i, bottom_edge);
            U.rho[ghost_B]  = U.rho[valid_B];
            U.rhou[ghost_B] = U.rhou[valid_B];
            U.rhov[ghost_B] = U.rhov[valid_B];
            U.E[ghost_B]    = U.E[valid_B];

            // Top boundary: copy the nearest interior cell
            int ghost_T = U.index(i, top_edge + 1 + j);
            int valid_T = U.index(i, top_edge);
            U.rho[ghost_T]  = U.rho[valid_T];
            U.rhou[ghost_T] = U.rhou[valid_T];
            U.rhov[ghost_T] = U.rhov[valid_T];
            U.E[ghost_T]    = U.E[valid_T];
        }
    }
}