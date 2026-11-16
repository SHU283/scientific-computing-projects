#pragma once

#include <string>
#include <string_view>

namespace fp {

// Runtime configuration object; each member stores one runtime setting.
struct RuntimeConfig {
    // Default backend is "serial", meaning the serial CPU path.
    std::string backend = "serial";
    // Study mode used for this run.
    std::string study = "single";
    // Numerical method selected for the run.
    std::string method = "godunov";
    // Test problem selected for the run.
    std::string problem = "sod";
    // Boundary-condition type.
    std::string boundary = "outflow";
    // Default grid size is 200 cells in x and 1 cell in y.
    int nx = 200;
    int ny = 1;
    // Maximum number of time steps allowed.
    int max_steps = 20000;
    // Computational domain range in the x direction is [0, 1].
    double x_min = 0.0;
    double x_max = 1.0;
    // Default CFL number is 0.4.
    double cfl = 0.4;
    // Default final time is 0.2.
    double final_time = 0.2;
    // Default ideal-gas heat-capacity ratio gamma is 1.4.
    double gamma = 1.4;
    // Default velocity for advection tests is 0.73.
    double advection_velocity = 0.73;
    // Path where the result file is written.
    std::string output_path = "results/stageA/sod_1d_godunov.csv";
};

// Create a default runtime configuration.
RuntimeConfig make_default_config();
// Parse the runtime configuration from command-line arguments.
RuntimeConfig parse_runtime_config(int argc, char** argv);
// Print the current runtime configuration.
void print_config(const RuntimeConfig& config);
// Check whether a problem name is supported.
bool is_supported_problem(std::string_view problem);
// Check whether a method name is supported.
bool is_supported_method(std::string_view method);
// Generate the Stage A output path from the current configuration.
std::string stage_a_output_path(const RuntimeConfig& config);

}  // namespace fp

