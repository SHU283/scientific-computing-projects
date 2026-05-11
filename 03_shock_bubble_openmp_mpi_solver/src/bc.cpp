#include "bc.h"

// Copy one cell from isrc to idst
static inline void copy_cell(Field& U, int idst, int isrc) {
    U.rho[idst]  = U.rho[isrc];
    U.rhou[idst] = U.rhou[isrc];
    U.rhov[idst] = U.rhov[isrc];
    U.E[idst]    = U.E[isrc];
}

// Fill left and right ghost cells
void apply_bc_lr_transmissive(Field& U)
{
    int ng = U.ng;                                                                              
    int il = ng;                                                                                
    int ir = ng + U.nx - 1;

    // for all rows including ghosts (0..height-1), fill left/right ghost columns
    for (int j = 0; j < U.height; ++j) {
        for (int g = 1; g <= ng; ++g) {
            // left ghost: i = il - g copies i = il
            int id_dst = U.idx(il - g, j);                                                      
            int id_src = U.idx(il, j);
            copy_cell(U, id_dst, id_src);

            // right ghost: i = ir + g copies i = ir
            id_dst = U.idx(ir + g, j);
            id_src = U.idx(ir, j);
            copy_cell(U, id_dst, id_src);
        }
    }
}

// Fill top and bottom ghost cells
void apply_bc_tb_transmissive(Field& U)
{
    int ng = U.ng;
    int jb = ng;
    int jt = ng + U.ny - 1;

    // for all columns including ghosts (0..pitch-1), fill bottom/top ghost rows
    for (int i = 0; i < U.pitch; ++i) {
        for (int g = 1; g <= ng; ++g) {
            // bottom ghost: j = jb - g copies j = jb
            int id_dst = U.idx(i, jb - g);
            int id_src = U.idx(i, jb);
            copy_cell(U, id_dst, id_src);

            // top ghost: j = jt + g copies j = jt
            id_dst = U.idx(i, jt + g);
            id_src = U.idx(i, jt);
            copy_cell(U, id_dst, id_src);
        }
    }
}

// Combine left/right and top/bottom boundary updates
void apply_bc_transmissive(Field& U)
{
    apply_bc_lr_transmissive(U);
    apply_bc_tb_transmissive(U);
} 
