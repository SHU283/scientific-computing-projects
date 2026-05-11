#pragma once
#include <vector> 

// Storage for a 2D field with ghost cells
struct Field {
    // Number of interior cells in x/y, ghost-cell count,
    // and total storage dimensions including ghosts
    int nx, ny;          
    int ng;             
    int pitch;          
    int height;        

    // Conserved variable arrays
    std::vector<double> rho, rhou, rhov, E;

    Field(int nx_, int ny_, int ng_=2){
        nx=nx_;
        ny=ny_;
        ng=ng_;
        pitch  = nx + 2*ng;       
        height = ny + 2*ng;        
        int N = pitch * height;   
        rho.resize(N); rhou.resize(N); rhov.resize(N); E.resize(N);   
    }

    inline int index(int i, int j) const {
        return i + pitch * j;
    }
};
