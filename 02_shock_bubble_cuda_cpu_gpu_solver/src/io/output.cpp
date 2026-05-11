#include "core/field.h"
#include <fstream>
#include <string>
#include "physics/state_convert.h"

// Create file to store density data, skip ghost
void write_rho_csv(const Field& U, const std::string& path)           
{
    std::ofstream f(path);                             
    for (int j=U.ng; j<U.ng+U.ny; ++j) {             
        for (int i=U.ng; i<U.ng+U.nx; ++i) { 
            f << U.rho[U.index(i,j)];                 
            if (i != U.ng + U.nx - 1) f << ",";     
        }
        f << "\n";
    }
}

// Create file to store pressure data, skip ghost
void write_p_csv(const Field& U, const std::string& path, double gamma) 
{
    std::ofstream f(path);
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            int id = U.index(i, j);
            double p = E_to_p(
                U.rho[id],
                U.rhou[id],
                U.rhov[id],
                U.E[id],
                gamma
            );
            f << p;
            if (i != U.ng + U.nx - 1) f << ",";
        }
        f << "\n";
    }
}
