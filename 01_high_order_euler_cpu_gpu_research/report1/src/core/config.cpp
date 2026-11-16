#include "config.h"

#include <iostream>
#include <stdexcept>

namespace {

bool starts_with(const std::string& text, const std::string& prefix) {
    return text.rfind(prefix, 0) == 0;
}

}  // namespace

namespace fp {

namespace {

// Stage A outputs are kept under a single result directory.
std::string stage_a_base_dir() {
    return "results/stageA/";
}

std::string make_stage_a_output_path(std::string_view problem, std::string_view method) {
    // 2D studies use the problem name directly because nx and ny both matter.
    if (problem == "gaussian_advection_2d" || problem == "quadrant_riemann_2d" || problem == "shock_bubble_2d") {
        return stage_a_base_dir() + std::string(problem) + "_" + std::string(method) + ".csv";
    }
    // 1D outputs are labelled explicitly to avoid confusion with 2D variants.
    return stage_a_base_dir() + std::string(problem) + "_1d_" + std::string(method) + ".csv";
}

// Apply a single --key=value command-line argument to the runtime configuration.
void apply_argument(RuntimeConfig& config, const std::string& arg) {
    if (!starts_with(arg, "--")) {
        throw std::runtime_error("Unsupported argument format: " + arg);
    }

    const std::size_t pos = arg.find('=');
    if (pos == std::string::npos) {
        throw std::runtime_error("Argument must use --key=value format: " + arg);
    }

    const std::string key = arg.substr(2, pos - 2);
    const std::string value = arg.substr(pos + 1);

    if (key == "backend") {
        config.backend = value;
    } else if (key == "study") {
        config.study = value;
    } else if (key == "boundary") {
        config.boundary = value;
    } else if (key == "method") {
        config.method = value;
        config.output_path = make_stage_a_output_path(config.problem, config.method);
    } else if (key == "problem") {
        config.problem = value;
        if (value == "smooth_advection" || value == "gaussian_advection" || value == "gaussian_advection_2d") {
            // Smooth advection benchmarks are periodic and use a longer final time.
            config.boundary = "periodic";
            config.final_time = 0.3;
            // If the user selects the 2D Gaussian problem without ny, mirror nx.
            if (value == "gaussian_advection_2d" && config.ny == 1) {
                config.ny = config.nx;       
            } else if (value == "sod") {                
                config.boundary = "outflow";
                config.final_time = 0.2;
            }
        }
        config.output_path = make_stage_a_output_path(config.problem, config.method);
    } else if (key == "nx") {
        config.nx = std::stoi(value);
    } else if (key == "ny") {
        config.ny = std::stoi(value);
    } else if (key == "max_steps") {
        config.max_steps = std::stoi(value);
    } else if (key == "x_min") {
        config.x_min = std::stod(value);
    } else if (key == "x_max") {
        config.x_max = std::stod(value);
    } else if (key == "cfl") {
        config.cfl = std::stod(value);
    } else if (key == "final_time") {
        config.final_time = std::stod(value);
    } else if (key == "gamma") {
        config.gamma = std::stod(value);
    } else if (key == "advection_velocity") {
        config.advection_velocity = std::stod(value);
    } else if (key == "output") {
        config.output_path = value;
    } else {
        throw std::runtime_error("Unknown argument key: " + key);
    }
}

}  // namespace

// Default runtime configuration used when no command-line overrides are given.
RuntimeConfig make_default_config() {
    RuntimeConfig config;
    config.output_path = stage_a_output_path(config);
    return config;
}

// Parse command-line arguments and validate the resulting configuration.
RuntimeConfig parse_runtime_config(int argc, char** argv) {
    RuntimeConfig config = make_default_config();

    for (int i = 1; i < argc; ++i) {
        apply_argument(config, argv[i]);
    }

    if (!is_supported_problem(config.problem)) {
        throw std::runtime_error("Unsupported problem: " + config.problem);
    }

    if (!is_supported_method(config.method)) {
        throw std::runtime_error("Unsupported method: " + config.method);
    }

    if (config.nx <= 0) {
        throw std::runtime_error("nx must be positive.");
    }

    if (config.x_max <= config.x_min) {
        throw std::runtime_error("x_max must be greater than x_min.");
    }

    return config;
}

// Print the resolved configuration before running a solver.
void print_config(const RuntimeConfig& config) {
    std::cout << "========== Final Project Solver ==========\n";
    std::cout << "backend    : " << config.backend << "\n";
    std::cout << "method     : " << config.method << "\n";
    std::cout << "problem    : " << config.problem << "\n";
    std::cout << "boundary   : " << config.boundary << "\n";
    std::cout << "grid       : " << config.nx << " x " << config.ny << "\n";
    std::cout << "cfl        : " << config.cfl << "\n";
    std::cout << "final_time : " << config.final_time << "\n";
    std::cout << "adv_vel    : " << config.advection_velocity << "\n";
    std::cout << "output     : " << config.output_path << "\n";
    std::cout << "==========================================\n";
}

bool is_supported_problem(std::string_view problem) {
    return problem == "sod" || problem == "smooth_advection" || problem == "gaussian_advection" ||
           problem == "gaussian_advection_2d" || problem == "quadrant_riemann_2d" ||
           problem == "shock_bubble_2d";
}

bool is_supported_method(std::string_view method) {
    return method == "godunov" || method == "muscl" || method == "weno2" || method == "weno3" ||
           method == "ader3";
}

// Default Stage A output path for the selected problem and method.
std::string stage_a_output_path(const RuntimeConfig& config) {
    return make_stage_a_output_path(config.problem, config.method);
}

}  // namespace fp

