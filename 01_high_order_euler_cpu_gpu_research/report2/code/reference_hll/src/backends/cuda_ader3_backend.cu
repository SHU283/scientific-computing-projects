#include "cuda_backend.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kDerivativeStep = 1.0e-3;
constexpr double kDirectionalDifferenceScale = 1.0e-6;
constexpr double kTroubledJumpThreshold = 0.15;
constexpr int kMaximumStepRetries = 8;

struct Ader3Polynomial2D {
    fp::Conservative2D constant;
    fp::Conservative2D slope_x;
    fp::Conservative2D slope_y;
    fp::Conservative2D curvature_x;
    fp::Conservative2D curvature_xy;
    fp::Conservative2D curvature_y;
};

__device__ fp::Conservative2D a3_add(
    const fp::Conservative2D& a,
    const fp::Conservative2D& b) {
    return {
        a.rho + b.rho,
        a.mom_x + b.mom_x,
        a.mom_y + b.mom_y,
        a.energy + b.energy,
    };
}

__device__ fp::Conservative2D a3_sub(
    const fp::Conservative2D& a,
    const fp::Conservative2D& b) {
    return {
        a.rho - b.rho,
        a.mom_x - b.mom_x,
        a.mom_y - b.mom_y,
        a.energy - b.energy,
    };
}

__device__ fp::Conservative2D a3_scale(
    double scalar,
    const fp::Conservative2D& value) {
    return {
        scalar * value.rho,
        scalar * value.mom_x,
        scalar * value.mom_y,
        scalar * value.energy,
    };
}

__device__ int a3_sample_index(int index, int size, int periodic) {
    if (periodic != 0) {
        index %= size;
        return index < 0 ? index + size : index;
    }
    if (index < 0) {
        return 0;
    }
    return index >= size ? size - 1 : index;
}

__device__ fp::Conservative2D a3_sample(
    const fp::Conservative2D* cells,
    int nx,
    int ny,
    int i,
    int j,
    int periodic) {
    const int ii = a3_sample_index(i, nx, periodic);
    const int jj = a3_sample_index(j, ny, periodic);
    return cells[
        static_cast<std::size_t>(jj) * static_cast<std::size_t>(nx) +
        static_cast<std::size_t>(ii)];
}

__device__ fp::Primitive2D a3_primitive(
    const fp::Conservative2D& state,
    double gamma) {
    const double inv_rho = 1.0 / state.rho;
    const double u = state.mom_x * inv_rho;
    const double v = state.mom_y * inv_rho;
    const double kinetic =
        0.5 * (state.mom_x * state.mom_x + state.mom_y * state.mom_y) *
        inv_rho;
    return {
        state.rho,
        u,
        v,
        (gamma - 1.0) * (state.energy - kinetic),
    };
}

__device__ bool a3_admissible(
    const fp::Conservative2D& state,
    double gamma) {
    if (!(state.rho > 0.0) || !isfinite(state.rho) ||
        !isfinite(state.mom_x) || !isfinite(state.mom_y) ||
        !isfinite(state.energy)) {
        return false;
    }
    const double kinetic =
        0.5 * (state.mom_x * state.mom_x + state.mom_y * state.mom_y) /
        state.rho;
    return (gamma - 1.0) * (state.energy - kinetic) > 0.0;
}

__device__ fp::Conservative2D a3_flux_x(
    const fp::Conservative2D& state,
    double gamma) {
    const fp::Primitive2D primitive = a3_primitive(state, gamma);
    return {
        state.mom_x,
        state.mom_x * primitive.u + primitive.p,
        state.mom_y * primitive.u,
        (state.energy + primitive.p) * primitive.u,
    };
}

__device__ fp::Conservative2D a3_flux_y(
    const fp::Conservative2D& state,
    double gamma) {
    const fp::Primitive2D primitive = a3_primitive(state, gamma);
    return {
        state.mom_y,
        state.mom_x * primitive.v,
        state.mom_y * primitive.v + primitive.p,
        (state.energy + primitive.p) * primitive.v,
    };
}

__device__ double a3_sound_speed(
    const fp::Primitive2D& state,
    double gamma) {
    return sqrt(gamma * state.p / state.rho);
}

__device__ double a3_signal_x(
    const fp::Conservative2D& state,
    double gamma) {
    const fp::Primitive2D primitive = a3_primitive(state, gamma);
    return fabs(primitive.u) + a3_sound_speed(primitive, gamma);
}

__device__ double a3_signal_y(
    const fp::Conservative2D& state,
    double gamma) {
    const fp::Primitive2D primitive = a3_primitive(state, gamma);
    return fabs(primitive.v) + a3_sound_speed(primitive, gamma);
}

__device__ fp::Conservative2D a3_rusanov_x(
    const fp::Conservative2D& left,
    const fp::Conservative2D& right,
    double gamma) {
    const double alpha = fmax(a3_signal_x(left, gamma), a3_signal_x(right, gamma));
    return a3_sub(
        a3_scale(0.5, a3_add(a3_flux_x(left, gamma), a3_flux_x(right, gamma))),
        a3_scale(0.5 * alpha, a3_sub(right, left)));
}

__device__ fp::Conservative2D a3_rusanov_y(
    const fp::Conservative2D& bottom,
    const fp::Conservative2D& top,
    double gamma) {
    const double alpha =
        fmax(a3_signal_y(bottom, gamma), a3_signal_y(top, gamma));
    return a3_sub(
        a3_scale(0.5, a3_add(a3_flux_y(bottom, gamma), a3_flux_y(top, gamma))),
        a3_scale(0.5 * alpha, a3_sub(top, bottom)));
}

__device__ fp::Conservative2D a3_evaluate(
    const Ader3Polynomial2D& polynomial,
    double xi,
    double eta) {
    fp::Conservative2D value = polynomial.constant;
    value = a3_add(value, a3_scale(xi, polynomial.slope_x));
    value = a3_add(value, a3_scale(eta, polynomial.slope_y));
    value = a3_add(
        value, a3_scale(0.5 * xi * xi, polynomial.curvature_x));
    value = a3_add(
        value, a3_scale(xi * eta, polynomial.curvature_xy));
    value = a3_add(
        value, a3_scale(0.5 * eta * eta, polynomial.curvature_y));
    return value;
}

__device__ Ader3Polynomial2D a3_reconstruct(
    const fp::Conservative2D* cells,
    int nx,
    int ny,
    int i,
    int j,
    int periodic) {
    const fp::Conservative2D center =
        a3_sample(cells, nx, ny, i, j, periodic);
    const fp::Conservative2D left =
        a3_sample(cells, nx, ny, i - 1, j, periodic);
    const fp::Conservative2D right =
        a3_sample(cells, nx, ny, i + 1, j, periodic);
    const fp::Conservative2D bottom =
        a3_sample(cells, nx, ny, i, j - 1, periodic);
    const fp::Conservative2D top =
        a3_sample(cells, nx, ny, i, j + 1, periodic);
    const fp::Conservative2D curvature_x =
        a3_add(a3_sub(right, a3_scale(2.0, center)), left);
    const fp::Conservative2D curvature_y =
        a3_add(a3_sub(top, a3_scale(2.0, center)), bottom);
    const fp::Conservative2D curvature_xy = a3_scale(
        0.25,
        a3_add(
            a3_sub(
                a3_sample(cells, nx, ny, i + 1, j + 1, periodic),
                a3_sample(cells, nx, ny, i + 1, j - 1, periodic)),
            a3_sub(
                a3_sample(cells, nx, ny, i - 1, j - 1, periodic),
                a3_sample(cells, nx, ny, i - 1, j + 1, periodic))));
    return {
        a3_sub(center, a3_scale(1.0 / 24.0, a3_add(curvature_x, curvature_y))),
        a3_scale(0.5, a3_sub(right, left)),
        a3_scale(0.5, a3_sub(top, bottom)),
        curvature_x,
        curvature_xy,
        curvature_y,
    };
}

__device__ bool a3_polynomial_admissible(
    const Ader3Polynomial2D& polynomial,
    double gamma) {
    const double points[5] = {
        -0.5 - 2.0 * kDerivativeStep,
        -0.5,
        0.0,
        0.5,
        0.5 + 2.0 * kDerivativeStep,
    };
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            if (!a3_admissible(
                    a3_evaluate(polynomial, points[i], points[j]), gamma)) {
                return false;
            }
        }
    }
    return true;
}

__device__ double a3_component_scale(const fp::Conservative2D& state) {
    return fmax(
        1.0,
        fmax(
            fmax(fabs(state.rho), fabs(state.mom_x)),
            fmax(fabs(state.mom_y), fabs(state.energy))));
}

__device__ bool a3_jacobian_vector_product(
    const fp::Conservative2D& state,
    const fp::Conservative2D& direction,
    double gamma,
    bool x_direction,
    fp::Conservative2D* output) {
    double epsilon =
        kDirectionalDifferenceScale *
        a3_component_scale(state) / a3_component_scale(direction);
    fp::Conservative2D plus = a3_add(state, a3_scale(epsilon, direction));
    fp::Conservative2D minus = a3_sub(state, a3_scale(epsilon, direction));
    for (int attempt = 0;
         attempt < 8 &&
         (!a3_admissible(plus, gamma) || !a3_admissible(minus, gamma));
         ++attempt) {
        epsilon *= 0.25;
        plus = a3_add(state, a3_scale(epsilon, direction));
        minus = a3_sub(state, a3_scale(epsilon, direction));
    }
    if (!a3_admissible(plus, gamma) || !a3_admissible(minus, gamma)) {
        return false;
    }
    const fp::Conservative2D flux_plus =
        x_direction ? a3_flux_x(plus, gamma) : a3_flux_y(plus, gamma);
    const fp::Conservative2D flux_minus =
        x_direction ? a3_flux_x(minus, gamma) : a3_flux_y(minus, gamma);
    *output = a3_scale(1.0 / (2.0 * epsilon), a3_sub(flux_plus, flux_minus));
    return true;
}

__device__ fp::Conservative2D a3_spatial_time_derivative(
    const Ader3Polynomial2D& polynomial,
    double xi,
    double eta,
    double cell_dx,
    double cell_dy,
    double gamma) {
    const double h = kDerivativeStep;
    const fp::Conservative2D x_difference = a3_sub(
        a3_flux_x(a3_evaluate(polynomial, xi + h, eta), gamma),
        a3_flux_x(a3_evaluate(polynomial, xi - h, eta), gamma));
    const fp::Conservative2D y_difference = a3_sub(
        a3_flux_y(a3_evaluate(polynomial, xi, eta + h), gamma),
        a3_flux_y(a3_evaluate(polynomial, xi, eta - h), gamma));
    return a3_add(
        a3_scale(-1.0 / (2.0 * h * cell_dx), x_difference),
        a3_scale(-1.0 / (2.0 * h * cell_dy), y_difference));
}

__device__ bool a3_second_time_derivative(
    const Ader3Polynomial2D& polynomial,
    double xi,
    double eta,
    double cell_dx,
    double cell_dy,
    double gamma,
    fp::Conservative2D* output) {
    const double h = kDerivativeStep;
    const fp::Conservative2D state_x_plus =
        a3_evaluate(polynomial, xi + h, eta);
    const fp::Conservative2D state_x_minus =
        a3_evaluate(polynomial, xi - h, eta);
    const fp::Conservative2D ut_x_plus = a3_spatial_time_derivative(
        polynomial, xi + h, eta, cell_dx, cell_dy, gamma);
    const fp::Conservative2D ut_x_minus = a3_spatial_time_derivative(
        polynomial, xi - h, eta, cell_dx, cell_dy, gamma);
    fp::Conservative2D a_ut_plus{};
    fp::Conservative2D a_ut_minus{};
    if (!a3_jacobian_vector_product(
            state_x_plus, ut_x_plus, gamma, true, &a_ut_plus) ||
        !a3_jacobian_vector_product(
            state_x_minus, ut_x_minus, gamma, true, &a_ut_minus)) {
        return false;
    }

    const fp::Conservative2D state_y_plus =
        a3_evaluate(polynomial, xi, eta + h);
    const fp::Conservative2D state_y_minus =
        a3_evaluate(polynomial, xi, eta - h);
    const fp::Conservative2D ut_y_plus = a3_spatial_time_derivative(
        polynomial, xi, eta + h, cell_dx, cell_dy, gamma);
    const fp::Conservative2D ut_y_minus = a3_spatial_time_derivative(
        polynomial, xi, eta - h, cell_dx, cell_dy, gamma);
    fp::Conservative2D b_ut_plus{};
    fp::Conservative2D b_ut_minus{};
    if (!a3_jacobian_vector_product(
            state_y_plus, ut_y_plus, gamma, false, &b_ut_plus) ||
        !a3_jacobian_vector_product(
            state_y_minus, ut_y_minus, gamma, false, &b_ut_minus)) {
        return false;
    }

    *output = a3_add(
        a3_scale(
            -1.0 / (2.0 * h * cell_dx),
            a3_sub(a_ut_plus, a_ut_minus)),
        a3_scale(
            -1.0 / (2.0 * h * cell_dy),
            a3_sub(b_ut_plus, b_ut_minus)));
    return true;
}

__device__ bool a3_predict(
    const Ader3Polynomial2D& polynomial,
    double xi,
    double eta,
    double tau,
    double cell_dx,
    double cell_dy,
    double gamma,
    fp::Conservative2D* output) {
    const fp::Conservative2D state = a3_evaluate(polynomial, xi, eta);
    const fp::Conservative2D ut = a3_spatial_time_derivative(
        polynomial, xi, eta, cell_dx, cell_dy, gamma);
    fp::Conservative2D utt{};
    if (!a3_second_time_derivative(
            polynomial, xi, eta, cell_dx, cell_dy, gamma, &utt)) {
        return false;
    }
    *output = a3_add(
        a3_add(state, a3_scale(tau, ut)),
        a3_scale(0.5 * tau * tau, utt));
    return a3_admissible(*output, gamma);
}

__device__ bool a3_interface_troubled(
    const fp::Conservative2D& first,
    const fp::Conservative2D& second,
    double gamma,
    bool normal_is_x) {
    const fp::Primitive2D left = a3_primitive(first, gamma);
    const fp::Primitive2D right = a3_primitive(second, gamma);
    const double rho_jump =
        fabs(right.rho - left.rho) / fmax(fmax(left.rho, right.rho), 1.0e-12);
    const double pressure_jump =
        fabs(right.p - left.p) / fmax(fmax(left.p, right.p), 1.0e-12);
    const double sound_scale =
        fmax(a3_sound_speed(left, gamma), a3_sound_speed(right, gamma));
    const double left_velocity = normal_is_x ? left.u : left.v;
    const double right_velocity = normal_is_x ? right.u : right.v;
    const double velocity_jump =
        fabs(right_velocity - left_velocity) / fmax(sound_scale, 1.0e-12);
    return fmax(fmax(rho_jump, pressure_jump), velocity_jump) >
           kTroubledJumpThreshold;
}

__device__ fp::Conservative2D a3_time_averaged_x_flux(
    const fp::Conservative2D* cells,
    int nx,
    int ny,
    int iface,
    int j,
    int periodic,
    double gamma,
    double dt,
    double cell_dx,
    double cell_dy) {
    constexpr double gauss = 0.28867513459481288225;
    const int left_i = iface - 1;
    const int right_i = iface;
    const fp::Conservative2D left_average =
        a3_sample(cells, nx, ny, left_i, j, periodic);
    const fp::Conservative2D right_average =
        a3_sample(cells, nx, ny, right_i, j, periodic);
    const Ader3Polynomial2D left_polynomial =
        a3_reconstruct(cells, nx, ny, left_i, j, periodic);
    const Ader3Polynomial2D right_polynomial =
        a3_reconstruct(cells, nx, ny, right_i, j, periodic);

    if (a3_interface_troubled(left_average, right_average, gamma, true) ||
        !a3_polynomial_admissible(left_polynomial, gamma) ||
        !a3_polynomial_admissible(right_polynomial, gamma)) {
        return a3_rusanov_x(left_average, right_average, gamma);
    }

    const double face_nodes[2] = {-gauss, gauss};
    const double time_nodes[2] = {
        dt * (0.5 - gauss),
        dt * (0.5 + gauss),
    };
    fp::Conservative2D average{};
    for (int face_node = 0; face_node < 2; ++face_node) {
        for (int time_node = 0; time_node < 2; ++time_node) {
            fp::Conservative2D left{};
            fp::Conservative2D right{};
            if (!a3_predict(
                    left_polynomial, 0.5, face_nodes[face_node],
                    time_nodes[time_node], cell_dx, cell_dy, gamma, &left) ||
                !a3_predict(
                    right_polynomial, -0.5, face_nodes[face_node],
                    time_nodes[time_node], cell_dx, cell_dy, gamma, &right)) {
                return a3_rusanov_x(left_average, right_average, gamma);
            }
            average = a3_add(
                average, a3_scale(0.25, a3_rusanov_x(left, right, gamma)));
        }
    }
    return average;
}

__device__ fp::Conservative2D a3_time_averaged_y_flux(
    const fp::Conservative2D* cells,
    int nx,
    int ny,
    int i,
    int jface,
    int periodic,
    double gamma,
    double dt,
    double cell_dx,
    double cell_dy) {
    constexpr double gauss = 0.28867513459481288225;
    const int bottom_j = jface - 1;
    const int top_j = jface;
    const fp::Conservative2D bottom_average =
        a3_sample(cells, nx, ny, i, bottom_j, periodic);
    const fp::Conservative2D top_average =
        a3_sample(cells, nx, ny, i, top_j, periodic);
    const Ader3Polynomial2D bottom_polynomial =
        a3_reconstruct(cells, nx, ny, i, bottom_j, periodic);
    const Ader3Polynomial2D top_polynomial =
        a3_reconstruct(cells, nx, ny, i, top_j, periodic);

    if (a3_interface_troubled(bottom_average, top_average, gamma, false) ||
        !a3_polynomial_admissible(bottom_polynomial, gamma) ||
        !a3_polynomial_admissible(top_polynomial, gamma)) {
        return a3_rusanov_y(bottom_average, top_average, gamma);
    }

    const double face_nodes[2] = {-gauss, gauss};
    const double time_nodes[2] = {
        dt * (0.5 - gauss),
        dt * (0.5 + gauss),
    };
    fp::Conservative2D average{};
    for (int face_node = 0; face_node < 2; ++face_node) {
        for (int time_node = 0; time_node < 2; ++time_node) {
            fp::Conservative2D bottom{};
            fp::Conservative2D top{};
            if (!a3_predict(
                    bottom_polynomial, face_nodes[face_node], 0.5,
                    time_nodes[time_node], cell_dx, cell_dy, gamma, &bottom) ||
                !a3_predict(
                    top_polynomial, face_nodes[face_node], -0.5,
                    time_nodes[time_node], cell_dx, cell_dy, gamma, &top)) {
                return a3_rusanov_y(bottom_average, top_average, gamma);
            }
            average = a3_add(
                average, a3_scale(0.25, a3_rusanov_y(bottom, top, gamma)));
        }
    }
    return average;
}

__device__ double a3_atomic_max(double* address, double value) {
    auto* address_as_ull =
        reinterpret_cast<unsigned long long int*>(address);
    unsigned long long int old = *address_as_ull;
    while (__longlong_as_double(old) < value) {
        const unsigned long long int assumed = old;
        old = atomicCAS(
            address_as_ull, assumed, __double_as_longlong(value));
        if (old == assumed) {
            break;
        }
    }
    return __longlong_as_double(old);
}

__global__ void a3_signal_speed_kernel(
    const fp::Conservative2D* cells,
    int count,
    double gamma,
    double* max_x,
    double* max_y) {
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) {
        return;
    }
    a3_atomic_max(max_x, a3_signal_x(cells[index], gamma));
    a3_atomic_max(max_y, a3_signal_y(cells[index], gamma));
}

__global__ void a3_flux_x_kernel(
    const fp::Conservative2D* cells,
    fp::Conservative2D* fluxes,
    int nx,
    int ny,
    int periodic,
    double gamma,
    double dt,
    double cell_dx,
    double cell_dy) {
    const int iface = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (iface > nx || j >= ny) {
        return;
    }
    fluxes[
        static_cast<std::size_t>(j) * static_cast<std::size_t>(nx + 1) +
        static_cast<std::size_t>(iface)] =
        a3_time_averaged_x_flux(
            cells, nx, ny, iface, j, periodic, gamma, dt, cell_dx, cell_dy);
}

__global__ void a3_flux_y_kernel(
    const fp::Conservative2D* cells,
    fp::Conservative2D* fluxes,
    int nx,
    int ny,
    int periodic,
    double gamma,
    double dt,
    double cell_dx,
    double cell_dy) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int jface = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || jface > ny) {
        return;
    }
    fluxes[
        static_cast<std::size_t>(jface) * static_cast<std::size_t>(nx) +
        static_cast<std::size_t>(i)] =
        a3_time_averaged_y_flux(
            cells, nx, ny, i, jface, periodic, gamma, dt, cell_dx, cell_dy);
}

__global__ void a3_update_kernel(
    const fp::Conservative2D* cells,
    const fp::Conservative2D* flux_x,
    const fp::Conservative2D* flux_y,
    fp::Conservative2D* next,
    int nx,
    int ny,
    double dt_over_dx,
    double dt_over_dy,
    double gamma,
    int* invalid_cells) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) {
        return;
    }
    const std::size_t cell_index =
        static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
        static_cast<std::size_t>(i);
    const std::size_t x_left =
        static_cast<std::size_t>(j) * static_cast<std::size_t>(nx + 1) +
        static_cast<std::size_t>(i);
    const std::size_t x_right = x_left + 1;
    const std::size_t y_bottom =
        static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
        static_cast<std::size_t>(i);
    const std::size_t y_top =
        static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(nx) +
        static_cast<std::size_t>(i);
    fp::Conservative2D updated = cells[cell_index];
    updated = a3_sub(
        updated,
        a3_scale(dt_over_dx, a3_sub(flux_x[x_right], flux_x[x_left])));
    updated = a3_sub(
        updated,
        a3_scale(dt_over_dy, a3_sub(flux_y[y_top], flux_y[y_bottom])));
    next[cell_index] = updated;
    if (!a3_admissible(updated, gamma)) {
        atomicAdd(invalid_cells, 1);
    }
}

}  // namespace

namespace fp {

SimulationResult2D cuda_run_ader3_2d(
    const RuntimeConfig& config,
    const Grid2D& grid,
    const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("2D initial state size does not match grid.");
    }

    const std::size_t cell_count =
        static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    const std::size_t cell_bytes = sizeof(Conservative2D) * cell_count;
    const std::size_t flux_x_count =
        static_cast<std::size_t>(grid.nx + 1) *
        static_cast<std::size_t>(grid.ny);
    const std::size_t flux_y_count =
        static_cast<std::size_t>(grid.nx) *
        static_cast<std::size_t>(grid.ny + 1);
    const double cell_dx = dx(grid);
    const double cell_dy = dy(grid);
    const double inv_dx = 1.0 / cell_dx;
    const double inv_dy = 1.0 / cell_dy;
    const int periodic = config.boundary == "periodic" ? 1 : 0;

    Conservative2D* d_field = nullptr;
    Conservative2D* d_next = nullptr;
    Conservative2D* d_flux_x = nullptr;
    Conservative2D* d_flux_y = nullptr;
    double* d_max_x = nullptr;
    double* d_max_y = nullptr;
    int* d_invalid_cells = nullptr;

    auto free_device_memory = [&]() {
        if (d_field) cudaFree(d_field);
        if (d_next) cudaFree(d_next);
        if (d_flux_x) cudaFree(d_flux_x);
        if (d_flux_y) cudaFree(d_flux_y);
        if (d_max_x) cudaFree(d_max_x);
        if (d_max_y) cudaFree(d_max_y);
        if (d_invalid_cells) cudaFree(d_invalid_cells);
    };

    if (cudaMalloc(&d_field, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(
            &d_flux_x, sizeof(Conservative2D) * flux_x_count) != cudaSuccess ||
        cudaMalloc(
            &d_flux_y, sizeof(Conservative2D) * flux_y_count) != cudaSuccess ||
        cudaMalloc(&d_max_x, sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_max_y, sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_invalid_cells, sizeof(int)) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMalloc failed in cuda_run_ader3_2d.");
    }

    if (cudaMemcpy(
            d_field, initial_state.cells.data(), cell_bytes,
            cudaMemcpyHostToDevice) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_ader3_2d.");
    }

    const dim3 block2d(16, 16);
    const dim3 cell_grid(
        (grid.nx + block2d.x - 1) / block2d.x,
        (grid.ny + block2d.y - 1) / block2d.y);
    const dim3 flux_x_grid(
        (grid.nx + 1 + block2d.x - 1) / block2d.x,
        (grid.ny + block2d.y - 1) / block2d.y);
    const dim3 flux_y_grid(
        (grid.nx + block2d.x - 1) / block2d.x,
        (grid.ny + 1 + block2d.y - 1) / block2d.y);
    const int block1d = 256;
    const int grid1d =
        static_cast<int>((cell_count + block1d - 1) / block1d);

    double time = 0.0;
    int step = 0;
    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemset(d_max_x, 0, sizeof(double)) != cudaSuccess ||
            cudaMemset(d_max_y, 0, sizeof(double)) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error(
                "CFL reduction reset failed in cuda_run_ader3_2d.");
        }
        a3_signal_speed_kernel<<<grid1d, block1d>>>(
            d_field, static_cast<int>(cell_count), config.gamma,
            d_max_x, d_max_y);

        double max_x = 0.0;
        double max_y = 0.0;
        if (cudaMemcpy(
                &max_x, d_max_x, sizeof(double),
                cudaMemcpyDeviceToHost) != cudaSuccess ||
            cudaMemcpy(
                &max_y, d_max_y, sizeof(double),
                cudaMemcpyDeviceToHost) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error(
                "CFL reduction copy failed in cuda_run_ader3_2d.");
        }
        double dt = config.cfl /
            (fmax(max_x, 1.0e-12) * inv_dx +
             fmax(max_y, 1.0e-12) * inv_dy);
        dt = std::min(dt, config.final_time - time);

        bool accepted = false;
        for (int retry = 0; retry <= kMaximumStepRetries; ++retry) {
            a3_flux_x_kernel<<<flux_x_grid, block2d>>>(
                d_field, d_flux_x, grid.nx, grid.ny, periodic,
                config.gamma, dt, cell_dx, cell_dy);
            a3_flux_y_kernel<<<flux_y_grid, block2d>>>(
                d_field, d_flux_y, grid.nx, grid.ny, periodic,
                config.gamma, dt, cell_dx, cell_dy);
            if (cudaMemset(d_invalid_cells, 0, sizeof(int)) != cudaSuccess) {
                free_device_memory();
                throw std::runtime_error(
                    "Validity reset failed in cuda_run_ader3_2d.");
            }
            a3_update_kernel<<<cell_grid, block2d>>>(
                d_field, d_flux_x, d_flux_y, d_next, grid.nx, grid.ny,
                dt * inv_dx, dt * inv_dy, config.gamma, d_invalid_cells);
            if (cudaDeviceSynchronize() != cudaSuccess) {
                free_device_memory();
                throw std::runtime_error(
                    "ADER-3 kernels failed in cuda_run_ader3_2d.");
            }

            int invalid_cells = 0;
            if (cudaMemcpy(
                    &invalid_cells, d_invalid_cells, sizeof(int),
                    cudaMemcpyDeviceToHost) != cudaSuccess) {
                free_device_memory();
                throw std::runtime_error(
                    "Validity copy failed in cuda_run_ader3_2d.");
            }
            if (invalid_cells == 0) {
                accepted = true;
                break;
            }
            dt *= 0.5;
        }
        if (!accepted) {
            free_device_memory();
            throw std::runtime_error(
                "CUDA ADER-3 could not produce an admissible update.");
        }

        std::swap(d_field, d_next);
        time += dt;
        ++step;
    }

    if (step >= config.max_steps && time < config.final_time) {
        free_device_memory();
        throw std::runtime_error(
            "CUDA ADER-3 reached max_steps before final_time.");
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field.nx = grid.nx;
    result.field.ny = grid.ny;
    result.field.cells.resize(cell_count);
    if (cudaMemcpy(
            result.field.cells.data(), d_field, cell_bytes,
            cudaMemcpyDeviceToHost) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error(
            "Final cudaMemcpy D2H failed in cuda_run_ader3_2d.");
    }
    free_device_memory();
    result.final_time = time;
    result.steps = step;
    return result;
}

}  // namespace fp
