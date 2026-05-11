#pragma once
#include <vector>
#include <cstddef>

struct Field {
    // The number of grid points in the x and y directions and the number of ghost points
    int nx, ny, ng;             
    // Total storage width and height
    int pitch, height;                           
    // Four conserved quantity fields                                                   
    std::vector<double> rho, rhou, rhov, E;                                                         

    // Constructor: set dimensions and allocate storage
    Field(int nx_, int ny_, int ng_);           
    // Map 2D indices to a 1D storage index
    inline int idx(int i, int j) const { return i + pitch * j; }                                  
    // Return the total number of grid cells, preventing overflow  
    inline std::size_t size() const { return static_cast<std::size_t>(pitch) * height; }            
};
