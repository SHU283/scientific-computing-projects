#pragma once

#include <string>
#include <string_view>

namespace fp {

struct RuntimeConfig {
    std::string backend = "serial";
    std::string study = "single";
    std::string method = "godunov";
    std::string problem = "sod";
    std::string boundary = "outflow";
    int nx = 200;
    int ny = 1;
    int max_steps = 20000;
    double x_min = 0.0;
    double x_max = 1.0;
    double cfl = 0.4;
    double final_time = 0.2;
    double gamma = 1.4;
    double advection_velocity = 0.73;
    std::string output_path = "results/stageA/sod_1d_godunov.csv";
};

RuntimeConfig make_default_config();
RuntimeConfig parse_runtime_config(int argc, char** argv);
void print_config(const RuntimeConfig& config);
bool is_supported_problem(std::string_view problem);
bool is_supported_method(std::string_view method);
std::string stage_a_output_path(const RuntimeConfig& config);

}  // namespace fp
