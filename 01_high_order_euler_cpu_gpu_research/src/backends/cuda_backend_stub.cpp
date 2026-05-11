#include "cuda_backend.h"

#include <stdexcept>

namespace fp {

bool cuda_backend_compiled() {
    return false;
}

bool cuda_backend_probe() {
    return false;
}

std::string cuda_backend_summary() {
    return "CUDA backend not compiled into this build.";
}

SimulationResult1D cuda_copy_1d_state(const RuntimeConfig&, const std::vector<Conservative1D>&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult1D cuda_run_godunov_1d(const RuntimeConfig&, const std::vector<Conservative1D>&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult1D cuda_run_muscl_1d(const RuntimeConfig&, const std::vector<Conservative1D>&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult1D cuda_run_weno2_1d(const RuntimeConfig&, const std::vector<Conservative1D>&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult2D cuda_run_godunov_2d(const RuntimeConfig&, const Grid2D&, const Field2D&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult2D cuda_run_muscl_2d(const RuntimeConfig&, const Grid2D&, const Field2D&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

SimulationResult2D cuda_run_weno2_2d(const RuntimeConfig&, const Grid2D&, const Field2D&) {
    throw std::runtime_error("CUDA backend not compiled into this build.");
}

}  // namespace fp
