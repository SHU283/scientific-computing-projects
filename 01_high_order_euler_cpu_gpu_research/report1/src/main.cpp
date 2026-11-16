#include "analysis.h"
#include "analysis2d.h"
#include "config.h"
#include "cuda_backend.h"
#include "csv_writer.h"
#include "csv_writer2d.h"
#include "problems.h"
#include "problems2d.h"
#include "solver_factory.h"
#include "godunov_solver_2d.h"
#include "muscl_solver_2d.h"
#include "weno2_solver_2d.h"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

// Centralized location for Stage A validation and convergence outputs.
std::string stage_a_dir() {
    return "results/stageA/";
}

// Build the 2D grid metadata from the runtime configuration.
fp::Grid2D make_grid_2d(const fp::RuntimeConfig& config) {
    fp::Grid2D grid;
    grid.nx = config.nx;
    grid.ny = config.ny;
    grid.x_min = config.x_min;
    grid.x_max = config.x_max;
    grid.y_min = config.x_min;
    grid.y_max = config.x_max;
    return grid;
}

int run_single_case(const fp::RuntimeConfig& config) {
    if (config.backend == "cuda") {
        fp::print_config(config);
        std::cout << fp::cuda_backend_summary() << "\n";

        if (!fp::cuda_backend_compiled()) {
            std::cout << "CUDA solver kernels are not enabled in this build yet.\n";
            return 0;
        }

        const bool ok = fp::cuda_backend_probe();
        std::cout << "cuda smoke test : " << (ok ? "passed" : "failed") << "\n";
        if (!ok) {
            return 1;
        }

        if ((config.method == "godunov" || config.method == "muscl" || config.method == "weno2") &&
            (config.problem == "sod" || config.problem == "smooth_advection" || config.problem == "gaussian_advection")) {
            const auto initial_state = fp::make_initial_condition(config);
            fp::SimulationResult1D result;
            if (config.method == "godunov") {
                result = fp::cuda_run_godunov_1d(config, initial_state);
            } else if (config.method == "muscl") {
                result = fp::cuda_run_muscl_1d(config, initial_state);
            } else {
                result = fp::cuda_run_weno2_1d(config, initial_state);
            }
            fp::write_result_csv(config.output_path, result, config.gamma);
            std::cout << "cuda 1D " << config.method << " solver : passed\n";
            std::cout << "time steps : " << result.steps << "\n";
            std::cout << "final time : " << result.final_time << "\n";
            std::cout << "output csv : " << config.output_path << "\n";
            if (config.problem == "smooth_advection" || config.problem == "gaussian_advection") {
                const auto exact = fp::make_exact_solution(config, result.final_time);
                const auto norms = fp::compute_density_error(result, exact);
                fp::print_error_norms(norms);
            }
            std::cout << "Note: this is currently a minimal CUDA 1D path with GPU reconstruction/flux/update kernels.\n";
            return 0;
        }

        if ((config.method == "godunov" || config.method == "muscl" || config.method == "weno2") &&
            config.problem == "gaussian_advection_2d") {
            // Build the 2D grid metadata from the runtime configuration.
            const fp::Grid2D grid = make_grid_2d(config);
            const auto initial_state = fp::make_gaussian_advection_2d_initial_condition(config, grid);
            fp::SimulationResult2D result;
            if (config.method == "godunov") {
                result = fp::cuda_run_godunov_2d(config, grid, initial_state);
            } else if (config.method == "muscl") {
                result = fp::cuda_run_muscl_2d(config, grid, initial_state);
            } else {
                result = fp::cuda_run_weno2_2d(config, grid, initial_state);
            }
            fp::write_result_csv_2d(config.output_path, result, config.gamma);

            // Print run information or diagnostics for the user.
            std::cout << "cuda 2D " << config.method << " solver : passed\n";
            std::cout << "time steps : " << result.steps << "\n";
            std::cout << "final time : " << result.final_time << "\n";
            std::cout << "output csv : " << config.output_path << "\n";

            const auto exact = fp::make_gaussian_advection_2d_exact_solution(config, grid, result.final_time);
            const auto norms = fp::compute_density_error_2d(result.field, exact);
            fp::print_error_norms_2d(norms);
            std::cout << "Note: this is currently a minimal CUDA 2D path with GPU flux/reconstruction/update kernels.\n";
            return 0;
        }

        if ((config.method == "godunov" || config.method == "muscl" || config.method == "weno2") &&
            (config.problem == "quadrant_riemann_2d" || config.problem == "shock_bubble_2d")) {
            // Build the 2D grid metadata from the runtime configuration.
            const fp::Grid2D grid = make_grid_2d(config);
            const auto initial_state = (config.problem == "quadrant_riemann_2d")
                                           ? fp::make_quadrant_riemann_2d_initial_condition(config, grid)
                                           : fp::make_shock_bubble_2d_initial_condition(config, grid);
            fp::SimulationResult2D result;
            if (config.method == "godunov") {
                result = fp::cuda_run_godunov_2d(config, grid, initial_state);
            } else if (config.method == "muscl") {
                result = fp::cuda_run_muscl_2d(config, grid, initial_state);
            } else {
                result = fp::cuda_run_weno2_2d(config, grid, initial_state);
            }
            fp::write_result_csv_2d(config.output_path, result, config.gamma);

            std::cout << "cuda 2D " << config.method << " solver : passed\n";
            std::cout << "time steps : " << result.steps << "\n";
            std::cout << "final time : " << result.final_time << "\n";
            std::cout << "output csv : " << config.output_path << "\n";
            std::cout << "2D discontinuous CUDA case finished. Please inspect the exported field data.\n";
            std::cout << "Note: this is currently a minimal CUDA 2D path for discontinuous cases.\n";
            return 0;
        }

        std::cout << "Note: full CUDA solver path is not implemented yet; this is only the backend scaffold.\n";
        return 0;
    }

    if (config.problem == "gaussian_advection_2d" || config.problem == "quadrant_riemann_2d" ||
        config.problem == "shock_bubble_2d") {
        // Print run information or diagnostics for the user.
        fp::print_config(config);

        const fp::Grid2D grid = make_grid_2d(config);

        fp::Field2D initial_state;
        if (config.problem == "gaussian_advection_2d") {
            initial_state = fp::make_gaussian_advection_2d_initial_condition(config, grid);
        } else if (config.problem == "shock_bubble_2d") {
            initial_state = fp::make_shock_bubble_2d_initial_condition(config, grid);
        } else {
            initial_state = fp::make_quadrant_riemann_2d_initial_condition(config, grid);
        }

        fp::SimulationResult2D result;
        if (config.method == "godunov") {
            result = fp::run_godunov_2d(config, grid, initial_state);
        } else if (config.method == "muscl") {
            result = fp::run_muscl_2d(config, grid, initial_state);
        } else if (config.method == "weno2") {
            result = fp::run_weno2_2d(config, grid, initial_state);
        } else {
            throw std::runtime_error("2D single run currently supports method=godunov, method=muscl, or method=weno2.");
        }
        fp::write_result_csv_2d(config.output_path, result, config.gamma);

        // Print run information or diagnostics for the user.
        std::cout << "\n";
        std::cout << "2D calculation finished.\n";
        std::cout << "time steps : " << result.steps << "\n";
        std::cout << "final time : " << result.final_time << "\n";
        std::cout << "output csv : " << config.output_path << "\n";
        if (config.problem == "gaussian_advection_2d") {
            const auto exact = fp::make_gaussian_advection_2d_exact_solution(config, grid, result.final_time);
            const auto norms = fp::compute_density_error_2d(result.field, exact);
            fp::print_error_norms_2d(norms);
        } else {
            std::cout << "2D discontinuous case finished. Please inspect the exported field data.\n";
        }
        return 0;
    }

    // Print run information or diagnostics for the user.
    fp::print_config(config);

    const auto initial_state = fp::make_initial_condition(config);
    // Parse, validate, or apply runtime configuration settings.
    const auto result = fp::run_selected_solver(config, initial_state);
    fp::write_result_csv(config.output_path, result, config.gamma);

    // Print run information or diagnostics for the user.
    std::cout << "\n";
    std::cout << "simulation complete.\n";
    std::cout << "time steps : " << result.steps << "\n";
    std::cout << "final time : " << result.final_time << "\n";
    std::cout << "output csv : " << config.output_path << "\n";
    if (config.problem == "smooth_advection" || config.problem == "gaussian_advection") {
        const auto exact = fp::make_exact_solution(config, result.final_time);
        const auto norms = fp::compute_density_error(result, exact);
        fp::print_error_norms(norms);
    }
    std::cout << "Tip: use --method=godunov, --method=muscl, or --method=weno2.\n";
    return 0;
}

// Convergence-study bookkeeping and reporting step.
int run_convergence_study(const fp::RuntimeConfig& base_config) {
    // Convergence-study bookkeeping and reporting step.
    if (base_config.problem == "gaussian_advection_2d") {
        if (base_config.method != "godunov" && base_config.method != "muscl" && base_config.method != "weno2") {
            throw std::runtime_error("2D convergence study currently supports only method=godunov, method=muscl, or method=weno2.");
        }

        // Build the 2D grid metadata from the runtime configuration.
        const std::vector<int> resolutions{25, 50, 100};
        std::vector<fp::ErrorNorms2D> norms;
        norms.reserve(resolutions.size());

        // Print run information or diagnostics for the user.
        std::cout << "Running 2D convergence study for method=" << base_config.method << "\n";

        // Convergence-study bookkeeping and reporting step.
        // Build the 2D grid metadata from the runtime configuration.
        for (const int n : resolutions) {
            // Parse, validate, or apply runtime configuration settings.
            fp::RuntimeConfig config = base_config;
            // Build the 2D grid metadata from the runtime configuration.
            config.nx = n;
            config.ny = n;
            config.output_path =
                stage_a_dir() + config.problem + "_" + config.method + "_n" + std::to_string(n) + ".csv";

            // Print run information or diagnostics for the user.
            std::cout << "\n";
            std::cout << "[case] n = " << n << " x " << n << "\n";

            // Build the 2D grid metadata from the runtime configuration.
            const fp::Grid2D grid = make_grid_2d(config);
            const auto initial_state = fp::make_gaussian_advection_2d_initial_condition(config, grid);
            fp::SimulationResult2D result;
            if (config.method == "godunov") {
                result = fp::run_godunov_2d(config, grid, initial_state);
            } else if (config.method == "muscl") {
                result = fp::run_muscl_2d(config, grid, initial_state);
            } else {
                result = fp::run_weno2_2d(config, grid, initial_state);
            }
            fp::write_result_csv_2d(config.output_path, result, config.gamma);

            const auto exact = fp::make_gaussian_advection_2d_exact_solution(config, grid, result.final_time);
            const auto error = fp::compute_density_error_2d(result.field, exact);
            norms.push_back(error);

            std::cout << "output csv : " << config.output_path << "\n";
            fp::print_error_norms_2d(error);
        }

        // Convergence-study bookkeeping and reporting step.
        // Convergence-study bookkeeping and reporting step.
        fp::print_convergence_table_2d(resolutions, norms);
        // Convergence-study bookkeeping and reporting step.
        const std::string summary_path =
            stage_a_dir() + base_config.problem + "_" + base_config.method + "_convergence.csv";
        // Convergence-study bookkeeping and reporting step.
        fp::write_convergence_csv_2d(summary_path, resolutions, norms);
        std::cout << "\nsummary csv : " << summary_path << "\n";
        return 0;
    }

    // One-dimensional problem setup or update step.
    if (base_config.problem != "smooth_advection" && base_config.problem != "gaussian_advection") {
        throw std::runtime_error(
            "Convergence study currently supports only problem=smooth_advection or problem=gaussian_advection.");
    }

    // Convergence-study bookkeeping and reporting step.
    std::vector<int> resolutions{50, 100, 200, 400};
    if (base_config.method == "weno3") {
        // Convergence-study bookkeeping and reporting step.
        resolutions.push_back(800);
    }
    std::vector<fp::ErrorNorms1D> norms;
    norms.reserve(resolutions.size());

    // Convergence-study bookkeeping and reporting step.
    std::cout << "Running convergence study for method=" << base_config.method << "\n";

    // Convergence-study bookkeeping and reporting step.
    for (const int nx : resolutions) {
        // Parse, validate, or apply runtime configuration settings.
        fp::RuntimeConfig config = base_config;
        config.nx = nx;
        config.output_path =
            stage_a_dir() + config.problem + "_1d_" + config.method + "_nx" + std::to_string(nx) + ".csv";

        // Print run information or diagnostics for the user.
        std::cout << "\n";
        std::cout << "[case] nx = " << nx << "\n";

        const auto initial_state = fp::make_initial_condition(config);
        const auto result = fp::run_selected_solver(config, initial_state);
        fp::write_result_csv(config.output_path, result, config.gamma);

        const auto exact = fp::make_exact_solution(config, result.final_time);
        const auto error = fp::compute_density_error(result, exact);
        norms.push_back(error);

        std::cout << "output csv : " << config.output_path << "\n";
        fp::print_error_norms(error);
    }

    // Convergence-study bookkeeping and reporting step.
    fp::print_convergence_table(resolutions, norms);
    const std::string summary_path =
        stage_a_dir() + base_config.problem + "_1d_" + base_config.method + "_convergence.csv";
    // Convergence-study bookkeeping and reporting step.
    fp::write_convergence_csv(summary_path, resolutions, norms);
    std::cout << "\nsummary csv : " << summary_path << "\n";
    return 0;
}

}  // namespace

// Parse, validate, or apply runtime configuration settings.
int main(int argc, char** argv) {
    try {
        // Parse, validate, or apply runtime configuration settings.
        // Parse, validate, or apply runtime configuration settings.
        const fp::RuntimeConfig config = fp::parse_runtime_config(argc, argv);

        if (config.study == "single") {
            return run_single_case(config);
        }

        // Convergence-study bookkeeping and reporting step.
        if (config.study == "convergence") {
            return run_convergence_study(config);
        }

        throw std::runtime_error("Unsupported study mode: " + config.study);
    } catch (const std::exception& ex) {
        // Print run information or diagnostics for the user.
        std::cerr << "program failed: " << ex.what() << "\n";
        std::cerr << "example: ./final_project --study=convergence --problem=smooth_advection --method=muscl\n";
        return 1;
    }
}


