#include "output.h"
#include "eos.h"
#include <fstream>

// Output the density field rho
void write_rho_csv(const Field& U, const std::string& path)
{
    std::ofstream f(path);
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            f << U.rho[U.idx(i,j)];
            if (i != U.ng + U.nx - 1) f << ",";
        }
        f << "\n";
    }
}

// Output the pressure field p
void write_p_csv(const Field& U, const Params& P, const std::string& path)
{
    std::ofstream f(path);
    for (int j = U.ng; j < U.ng + U.ny; ++j) {
        for (int i = U.ng; i < U.ng + U.nx; ++i) {
            int id = U.idx(i,j);
            double p = pressure_from_cons(U.rho[id], U.rhou[id], U.rhov[id], U.E[id], P.gamma);
            f << p;
            if (i != U.ng + U.nx - 1) f << ",";
        }
        f << "\n";
    }
}
