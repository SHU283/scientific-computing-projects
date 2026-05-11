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

std::string stage_a_base_dir() {
    return "results/stageA/";
}

bool is_2d_problem(std::string_view problem) {
    // 2D cases need ny-aware output names and route through the 2D solver path.
    return problem == "gaussian_advection_2d" || problem == "quadrant_riemann_2d" ||
           problem == "shock_bubble_2d";
}

bool is_periodic_problem(std::string_view problem) {
    // Smooth transport benchmarks use periodic boundaries so they have a clean
    // exact solution for error and convergence studies.
    return problem == "smooth_advection" || problem == "gaussian_advection" ||
           problem == "gaussian_advection_2d";
}

std::string make_stage_a_output_path(std::string_view problem, std::string_view method) {
    // Keep generated outputs grouped by benchmark and method for post-processing scripts.
    if (is_2d_problem(problem)) {
        return stage_a_base_dir() + std::string(problem) + "_" + std::string(method) + ".csv";
    }
    return stage_a_base_dir() + std::string(problem) + "_1d_" + std::string(method) + ".csv";
}

void apply_problem_defaults(RuntimeConfig& config) {
    // Periodic advection cases are accuracy benchmarks; discontinuous flow cases use outflow boundaries.
    if (is_periodic_problem(config.problem)) {
        config.boundary = "periodic";
        config.final_time = 0.3;
    } else {
        config.boundary = "outflow";
        config.final_time = 0.2;
    }

    if (config.problem == "gaussian_advection_2d" && config.ny == 1) {
        config.ny = config.nx;
    }
}

void apply_argument(RuntimeConfig& config, const std::string& arg) {
    // The command line intentionally uses only --key=value pairs to keep scripts simple and reproducible.
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
        apply_problem_defaults(config);
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

RuntimeConfig make_default_config() {
    RuntimeConfig config;
    config.output_path = stage_a_output_path(config);
    return config;
}

RuntimeConfig parse_runtime_config(int argc, char** argv) {
    // Build a complete validated configuration from defaults plus CLI overrides.
    // The extra user_set_ny flag preserves explicit rectangular-grid requests.
    RuntimeConfig config = make_default_config();
    bool user_set_ny = false;

    for (int i = 1; i < argc; ++i) {
        if (starts_with(argv[i], "--ny=")) {
            user_set_ny = true;
        }
        apply_argument(config, argv[i]);
    }

    if (!is_supported_problem(config.problem)) {
        throw std::runtime_error("Unsupported problem: " + config.problem);
    }

    if (!is_supported_method(config.method)) {
        throw std::runtime_error("Unsupported method: " + config.method);
    }

    if (config.problem == "gaussian_advection_2d" && !user_set_ny) {
        config.ny = config.nx;
    }

    if (config.nx <= 0 || config.ny <= 0) {
        throw std::runtime_error("nx and ny must be positive.");
    }

    if (config.x_max <= config.x_min) {
        throw std::runtime_error("x_max must be greater than x_min.");
    }

    return config;
}

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

std::string stage_a_output_path(const RuntimeConfig& config) {
    return make_stage_a_output_path(config.problem, config.method);
}

}  // namespace fp
