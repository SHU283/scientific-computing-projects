#include "csv_writer.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "euler1d.h"

namespace fp {

void write_result_csv(const std::string& path, const SimulationResult1D& result, double gamma) {
    if (path == "none" || path == "/dev/null") {
        return;
    }

    const std::filesystem::path output_path(path);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output csv file: " + path);
    }

    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "x,rho,u,p,energy\n";

    for (std::size_t i = 0; i < result.cells.size(); ++i) {
        const Primitive1D w = conservative_to_primitive(result.cells[i], gamma);
        out << result.x[i] << ','
            << result.cells[i].rho << ','
            << w.u << ','
            << w.p << ','
            << result.cells[i].energy << '\n';
    }
}

}  // namespace fp
