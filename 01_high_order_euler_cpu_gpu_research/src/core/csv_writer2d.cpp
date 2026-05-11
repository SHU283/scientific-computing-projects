#include "csv_writer2d.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace fp {

void write_result_csv_2d(const std::string& path, const SimulationResult2D& result, double gamma) {
    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open 2D output csv file: " + path);
    }

    out << "i,j,rho,u,v,p,energy\n";
    for (int j = 0; j < result.grid.ny; ++j) {
        for (int i = 0; i < result.grid.nx; ++i) {
            const Conservative2D& cell = result.field(i, j);
            const Primitive2D w = conservative_to_primitive(cell, gamma);
            out << i << ','
                << j << ','
                << cell.rho << ','
                << w.u << ','
                << w.v << ','
                << w.p << ','
                << cell.energy << '\n';
        }
    }
}

}  // namespace fp
