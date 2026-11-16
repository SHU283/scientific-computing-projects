#include "cuda_backend.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double kMinPositive = 1.0e-12;
constexpr double kWenoEps = 1.0e-6;
constexpr double kBackgroundPressure = 1.0;

__global__ void fp_cuda_smoke_kernel(int* out_flag) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        *out_flag = 1;
    }
}

__global__ void fp_cuda_copy_state_kernel(const fp::Conservative1D* in_state,
                                          fp::Conservative1D* out_state,
                                          int n) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        out_state[idx] = in_state[idx];
    }
}

__device__ fp::Primitive1D fp_cuda_conservative_to_primitive(const fp::Conservative1D& u, double gamma) {
    const double rho = fmax(u.rho, kMinPositive);
    const double velocity = u.mom / rho;
    const double kinetic = 0.5 * rho * velocity * velocity;
    const double pressure = fmax((gamma - 1.0) * (u.energy - kinetic), kMinPositive);
    return {rho, velocity, pressure};
}

__device__ double fp_cuda_sound_speed(const fp::Primitive1D& w, double gamma) {
    return sqrt(gamma * w.p / w.rho);
}

__device__ fp::Conservative1D fp_cuda_physical_flux(const fp::Conservative1D& u, double gamma) {
    const fp::Primitive1D w = fp_cuda_conservative_to_primitive(u, gamma);
    return {u.mom, u.mom * w.u + w.p, (u.energy + w.p) * w.u};
}

__device__ double fp_cuda_nonlinear_weight_z(double beta, double tau, double linear_weight) {
    const double ratio = tau / (kWenoEps + beta);
    return linear_weight * (1.0 + ratio * ratio);
}

__device__ double fp_cuda_mapped_weight(double omega, double linear_weight) {
    const double d = linear_weight;
    const double numerator = omega * (d + d * d - 3.0 * d * omega + omega * omega);
    const double denominator = d * d + omega * (1.0 - 2.0 * d);
    if (fabs(denominator) < 1.0e-14) {
        return omega;
    }
    return numerator / denominator;
}

__device__ double fp_cuda_weno2_left_component(double um1, double u0, double up1) {
    const double beta0 = (u0 - um1) * (u0 - um1);
    const double beta1 = (up1 - u0) * (up1 - u0);
    const double tau = fabs(beta0 - beta1);

    double alpha0 = fp_cuda_nonlinear_weight_z(beta0, tau, 0.5);
    double alpha1 = fp_cuda_nonlinear_weight_z(beta1, tau, 0.5);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = fp_cuda_mapped_weight(omega0, 0.5);
    omega1 = fp_cuda_mapped_weight(omega1, 0.5);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    const double q0 = -0.5 * um1 + 1.5 * u0;
    const double q1 = 0.5 * u0 + 0.5 * up1;
    return omega0 * q0 + omega1 * q1;
}

__device__ double fp_cuda_weno2_right_component(double u0, double up1, double up2) {
    const double beta0 = (up1 - u0) * (up1 - u0);
    const double beta1 = (up2 - up1) * (up2 - up1);
    const double tau = fabs(beta0 - beta1);

    double alpha0 = fp_cuda_nonlinear_weight_z(beta0, tau, 0.5);
    double alpha1 = fp_cuda_nonlinear_weight_z(beta1, tau, 0.5);
    double omega0 = alpha0 / (alpha0 + alpha1);
    double omega1 = alpha1 / (alpha0 + alpha1);
    omega0 = fp_cuda_mapped_weight(omega0, 0.5);
    omega1 = fp_cuda_mapped_weight(omega1, 0.5);
    const double omega_sum = omega0 + omega1;
    omega0 /= omega_sum;
    omega1 /= omega_sum;

    const double q0 = 0.5 * u0 + 0.5 * up1;
    const double q1 = 1.5 * up1 - 0.5 * up2;
    return omega0 * q0 + omega1 * q1;
}

__device__ double fp_cuda_weno3_left_component(double um1, double u0, double up1) {
    const double q0 = -0.5 * um1 + 1.5 * u0;
    const double q1 = 0.5 * u0 + 0.5 * up1;
    const double beta0 = (u0 - um1) * (u0 - um1);
    const double beta1 = (up1 - u0) * (up1 - u0);
    const double tau = fabs(beta0 - beta1);

    double omega0 = fp_cuda_nonlinear_weight_z(beta0, tau, 1.0 / 3.0);
    double omega1 = fp_cuda_nonlinear_weight_z(beta1, tau, 2.0 / 3.0);
    const double alpha_sum = omega0 + omega1;
    omega0 = fp_cuda_mapped_weight(omega0 / alpha_sum, 1.0 / 3.0);
    omega1 = fp_cuda_mapped_weight(omega1 / alpha_sum, 2.0 / 3.0);
    return (omega0 * q0 + omega1 * q1) / (omega0 + omega1);
}

__device__ double fp_cuda_weno3_right_component(double u0, double up1, double up2) {
    const double q0 = 0.5 * u0 + 0.5 * up1;
    const double q1 = 1.5 * up1 - 0.5 * up2;
    const double beta0 = (up1 - u0) * (up1 - u0);
    const double beta1 = (up2 - up1) * (up2 - up1);
    const double tau = fabs(beta0 - beta1);

    double omega0 = fp_cuda_nonlinear_weight_z(beta0, tau, 2.0 / 3.0);
    double omega1 = fp_cuda_nonlinear_weight_z(beta1, tau, 1.0 / 3.0);
    const double alpha_sum = omega0 + omega1;
    omega0 = fp_cuda_mapped_weight(omega0 / alpha_sum, 2.0 / 3.0);
    omega1 = fp_cuda_mapped_weight(omega1 / alpha_sum, 1.0 / 3.0);
    return (omega0 * q0 + omega1 * q1) / (omega0 + omega1);
}

__device__ fp::Conservative1D fp_cuda_hll_flux(const fp::Conservative1D& left,
                                               const fp::Conservative1D& right,
                                               double gamma) {
    const fp::Primitive1D wl = fp_cuda_conservative_to_primitive(left, gamma);
    const fp::Primitive1D wr = fp_cuda_conservative_to_primitive(right, gamma);
    const double cl = fp_cuda_sound_speed(wl, gamma);
    const double cr = fp_cuda_sound_speed(wr, gamma);

    const double s_l = fmin(wl.u - cl, wr.u - cr);
    const double s_r = fmax(wl.u + cl, wr.u + cr);

    const fp::Conservative1D f_l = fp_cuda_physical_flux(left, gamma);
    const fp::Conservative1D f_r = fp_cuda_physical_flux(right, gamma);

    if (s_l >= 0.0) {
        return f_l;
    }
    if (s_r <= 0.0) {
        return f_r;
    }

    const double denom = s_r - s_l;
    return {(s_r * f_l.rho - s_l * f_r.rho + s_l * s_r * (right.rho - left.rho)) / denom,
            (s_r * f_l.mom - s_l * f_r.mom + s_l * s_r * (right.mom - left.mom)) / denom,
            (s_r * f_l.energy - s_l * f_r.energy + s_l * s_r * (right.energy - left.energy)) / denom};
}

__device__ fp::Primitive2D fp_cuda_conservative_to_primitive_2d(const fp::Conservative2D& u, double gamma) {
    const double rho = fmax(u.rho, kMinPositive);
    const double vel_x = u.mom_x / rho;
    const double vel_y = u.mom_y / rho;
    const double kinetic = 0.5 * rho * (vel_x * vel_x + vel_y * vel_y);
    const double pressure = fmax((gamma - 1.0) * (u.energy - kinetic), kMinPositive);
    return {rho, vel_x, vel_y, pressure};
}

__device__ fp::Conservative2D fp_cuda_physical_flux_x_2d(const fp::Conservative2D& u, double gamma) {
    const fp::Primitive2D w = fp_cuda_conservative_to_primitive_2d(u, gamma);
    return {u.mom_x, u.mom_x * w.u + w.p, u.mom_y * w.u, (u.energy + w.p) * w.u};
}

__device__ fp::Conservative2D fp_cuda_physical_flux_y_2d(const fp::Conservative2D& u, double gamma) {
    const fp::Primitive2D w = fp_cuda_conservative_to_primitive_2d(u, gamma);
    return {u.mom_y, u.mom_x * w.v, u.mom_y * w.v + w.p, (u.energy + w.p) * w.v};
}

__device__ double fp_cuda_sound_speed_2d(const fp::Primitive2D& w, double gamma) {
    return sqrt(gamma * w.p / w.rho);
}

__device__ double fp_cuda_max_signal_speed_x_2d(const fp::Conservative2D& u, double gamma) {
    const fp::Primitive2D w = fp_cuda_conservative_to_primitive_2d(u, gamma);
    return fabs(w.u) + fp_cuda_sound_speed_2d(w, gamma);
}

__device__ double fp_cuda_max_signal_speed_y_2d(const fp::Conservative2D& u, double gamma) {
    const fp::Primitive2D w = fp_cuda_conservative_to_primitive_2d(u, gamma);
    return fabs(w.v) + fp_cuda_sound_speed_2d(w, gamma);
}

// F = 0.5*(F_L + F_R) - 0.5*alpha*(U_R - U_L).
__device__ fp::Conservative2D fp_cuda_rusanov_flux_x_2d(const fp::Conservative2D& left,
                                                        const fp::Conservative2D& right,
                                                        double gamma) {
    const fp::Conservative2D f_left = fp_cuda_physical_flux_x_2d(left, gamma);
    const fp::Conservative2D f_right = fp_cuda_physical_flux_x_2d(right, gamma);
    const double alpha = fmax(fp_cuda_max_signal_speed_x_2d(left, gamma), fp_cuda_max_signal_speed_x_2d(right, gamma));
    return {0.5 * (f_left.rho + f_right.rho) - 0.5 * alpha * (right.rho - left.rho),
            0.5 * (f_left.mom_x + f_right.mom_x) - 0.5 * alpha * (right.mom_x - left.mom_x),
            0.5 * (f_left.mom_y + f_right.mom_y) - 0.5 * alpha * (right.mom_y - left.mom_y),
            0.5 * (f_left.energy + f_right.energy) - 0.5 * alpha * (right.energy - left.energy)};
}

__device__ fp::Conservative2D fp_cuda_rusanov_flux_y_2d(const fp::Conservative2D& bottom,
                                                        const fp::Conservative2D& top,
                                                        double gamma) {
    const fp::Conservative2D f_bottom = fp_cuda_physical_flux_y_2d(bottom, gamma);
    const fp::Conservative2D f_top = fp_cuda_physical_flux_y_2d(top, gamma);
    const double alpha =
        fmax(fp_cuda_max_signal_speed_y_2d(bottom, gamma), fp_cuda_max_signal_speed_y_2d(top, gamma));
    return {0.5 * (f_bottom.rho + f_top.rho) - 0.5 * alpha * (top.rho - bottom.rho),
            0.5 * (f_bottom.mom_x + f_top.mom_x) - 0.5 * alpha * (top.mom_x - bottom.mom_x),
            0.5 * (f_bottom.mom_y + f_top.mom_y) - 0.5 * alpha * (top.mom_y - bottom.mom_y),
            0.5 * (f_bottom.energy + f_top.energy) - 0.5 * alpha * (top.energy - bottom.energy)};
}

__device__ double fp_cuda_minmod(double a, double b);

__device__ fp::Conservative2D fp_cuda_limited_slope_x_2d(const fp::Conservative2D& left,
                                                         const fp::Conservative2D& center,
                                                         const fp::Conservative2D& right) {
    return {fp_cuda_minmod(center.rho - left.rho, right.rho - center.rho),
            fp_cuda_minmod(center.mom_x - left.mom_x, right.mom_x - center.mom_x),
            fp_cuda_minmod(center.mom_y - left.mom_y, right.mom_y - center.mom_y),
            fp_cuda_minmod(center.energy - left.energy, right.energy - center.energy)};
}

__device__ fp::Conservative2D fp_cuda_limited_slope_y_2d(const fp::Conservative2D& bottom,
                                                         const fp::Conservative2D& center,
                                                         const fp::Conservative2D& top) {
    return {fp_cuda_minmod(center.rho - bottom.rho, top.rho - center.rho),
            fp_cuda_minmod(center.mom_x - bottom.mom_x, top.mom_x - center.mom_x),
            fp_cuda_minmod(center.mom_y - bottom.mom_y, top.mom_y - center.mom_y),
            fp_cuda_minmod(center.energy - bottom.energy, top.energy - center.energy)};
}

__device__ fp::Conservative2D fp_cuda_sample_cell_2d(const fp::Conservative2D* cells,
                                                     int nx,
                                                     int ny,
                                                     int i,
                                                     int j,
                                                     int periodic_boundary) {
    if (periodic_boundary != 0) {
        while (i < 0) i += nx;
        while (i >= nx) i -= nx;
        while (j < 0) j += ny;
        while (j >= ny) j -= ny;
    } else {
        if (i < 0) i = 0;
        if (i >= nx) i = nx - 1;
        if (j < 0) j = 0;
        if (j >= ny) j = ny - 1;
    }

    return cells[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i)];
}

__device__ double fp_cuda_minmod(double a, double b) {
    if (a * b <= 0.0) {
        return 0.0;
    }
    return (fabs(a) < fabs(b)) ? a : b;
}

__device__ fp::Conservative1D fp_cuda_limited_slope(const fp::Conservative1D& left_diff,
                                                    const fp::Conservative1D& right_diff) {
    return {fp_cuda_minmod(left_diff.rho, right_diff.rho),
            fp_cuda_minmod(left_diff.mom, right_diff.mom),
            fp_cuda_minmod(left_diff.energy, right_diff.energy)};
}

__device__ fp::Conservative1D fp_cuda_sample_cell_1d(const fp::Conservative1D* cells,
                                                     int nx,
                                                     int i,
                                                     int periodic_boundary) {
    if (periodic_boundary != 0) {
        while (i < 0) i += nx;
        while (i >= nx) i -= nx;
    } else {
        if (i < 0) i = 0;
        if (i >= nx) i = nx - 1;
    }
    return cells[i];
}

__device__ double fp_cuda_sample_scalar_1d(const double* cells, int nx, int i, int periodic_boundary) {
    if (periodic_boundary != 0) {
        while (i < 0) i += nx;
        while (i >= nx) i -= nx;
    } else {
        if (i < 0) i = 0;
        if (i >= nx) i = nx - 1;
    }
    return cells[i];
}

__device__ fp::Conservative1D fp_cuda_weno2_left_state_1d(const fp::Conservative1D& um1,
                                                          const fp::Conservative1D& u0,
                                                          const fp::Conservative1D& up1) {
    return {fp_cuda_weno2_left_component(um1.rho, u0.rho, up1.rho),
            fp_cuda_weno2_left_component(um1.mom, u0.mom, up1.mom),
            fp_cuda_weno2_left_component(um1.energy, u0.energy, up1.energy)};
}

__device__ fp::Conservative1D fp_cuda_weno2_right_state_1d(const fp::Conservative1D& u0,
                                                           const fp::Conservative1D& up1,
                                                           const fp::Conservative1D& up2) {
    return {fp_cuda_weno2_right_component(u0.rho, up1.rho, up2.rho),
            fp_cuda_weno2_right_component(u0.mom, up1.mom, up2.mom),
            fp_cuda_weno2_right_component(u0.energy, up1.energy, up2.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno2_left_state_x_2d(const fp::Conservative2D* cells,
                                                            int nx,
                                                            int ny,
                                                            int i,
                                                            int j,
                                                            int periodic_boundary) {
    const fp::Conservative2D um1 = fp_cuda_sample_cell_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    return {fp_cuda_weno2_left_component(um1.rho, u0.rho, up1.rho),
            fp_cuda_weno2_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
            fp_cuda_weno2_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
            fp_cuda_weno2_left_component(um1.energy, u0.energy, up1.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno2_right_state_x_2d(const fp::Conservative2D* cells,
                                                             int nx,
                                                             int ny,
                                                             int i,
                                                             int j,
                                                             int periodic_boundary) {
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    const fp::Conservative2D up2 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 2, j, periodic_boundary);
    return {fp_cuda_weno2_right_component(u0.rho, up1.rho, up2.rho),
            fp_cuda_weno2_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
            fp_cuda_weno2_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
            fp_cuda_weno2_right_component(u0.energy, up1.energy, up2.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno2_left_state_y_2d(const fp::Conservative2D* cells,
                                                            int nx,
                                                            int ny,
                                                            int i,
                                                            int j,
                                                            int periodic_boundary) {
    const fp::Conservative2D um1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    return {fp_cuda_weno2_left_component(um1.rho, u0.rho, up1.rho),
            fp_cuda_weno2_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
            fp_cuda_weno2_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
            fp_cuda_weno2_left_component(um1.energy, u0.energy, up1.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno2_right_state_y_2d(const fp::Conservative2D* cells,
                                                             int nx,
                                                             int ny,
                                                             int i,
                                                             int j,
                                                             int periodic_boundary) {
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    const fp::Conservative2D up2 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 2, periodic_boundary);
    return {fp_cuda_weno2_right_component(u0.rho, up1.rho, up2.rho),
            fp_cuda_weno2_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
            fp_cuda_weno2_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
             fp_cuda_weno2_right_component(u0.energy, up1.energy, up2.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno3_left_state_x_2d(const fp::Conservative2D* cells,
                                                            int nx,
                                                            int ny,
                                                            int i,
                                                            int j,
                                                            int periodic_boundary) {
    const fp::Conservative2D um1 = fp_cuda_sample_cell_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    return {fp_cuda_weno3_left_component(um1.rho, u0.rho, up1.rho),
            fp_cuda_weno3_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
            fp_cuda_weno3_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
            fp_cuda_weno3_left_component(um1.energy, u0.energy, up1.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno3_right_state_x_2d(const fp::Conservative2D* cells,
                                                             int nx,
                                                             int ny,
                                                             int i,
                                                             int j,
                                                             int periodic_boundary) {
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    const fp::Conservative2D up2 = fp_cuda_sample_cell_2d(cells, nx, ny, i + 2, j, periodic_boundary);
    return {fp_cuda_weno3_right_component(u0.rho, up1.rho, up2.rho),
            fp_cuda_weno3_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
            fp_cuda_weno3_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
            fp_cuda_weno3_right_component(u0.energy, up1.energy, up2.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno3_left_state_y_2d(const fp::Conservative2D* cells,
                                                            int nx,
                                                            int ny,
                                                            int i,
                                                            int j,
                                                            int periodic_boundary) {
    const fp::Conservative2D um1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    return {fp_cuda_weno3_left_component(um1.rho, u0.rho, up1.rho),
            fp_cuda_weno3_left_component(um1.mom_x, u0.mom_x, up1.mom_x),
            fp_cuda_weno3_left_component(um1.mom_y, u0.mom_y, up1.mom_y),
            fp_cuda_weno3_left_component(um1.energy, u0.energy, up1.energy)};
}

__device__ fp::Conservative2D fp_cuda_weno3_right_state_y_2d(const fp::Conservative2D* cells,
                                                             int nx,
                                                             int ny,
                                                             int i,
                                                             int j,
                                                             int periodic_boundary) {
    const fp::Conservative2D u0 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D up1 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    const fp::Conservative2D up2 = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 2, periodic_boundary);
    return {fp_cuda_weno3_right_component(u0.rho, up1.rho, up2.rho),
            fp_cuda_weno3_right_component(u0.mom_x, up1.mom_x, up2.mom_x),
            fp_cuda_weno3_right_component(u0.mom_y, up1.mom_y, up2.mom_y),
            fp_cuda_weno3_right_component(u0.energy, up1.energy, up2.energy)};
}

__device__ bool fp_cuda_is_admissible_2d(const fp::Conservative2D& state, double gamma) {
    if (!(state.rho > 0.0) || !isfinite(state.rho) || !isfinite(state.mom_x) ||
        !isfinite(state.mom_y) || !isfinite(state.energy)) {
        return false;
    }
    const double kinetic =
        0.5 * (state.mom_x * state.mom_x + state.mom_y * state.mom_y) / state.rho;
    return (gamma - 1.0) * (state.energy - kinetic) > 0.0;
}

__global__ void fp_cuda_godunov_flux_kernel(const fp::Conservative1D* cells,
                                            fp::Conservative1D* fluxes,
                                            int nx,
                                            double gamma,
                                            int periodic_boundary) {
    const int iface = blockIdx.x * blockDim.x + threadIdx.x;
    if (iface > nx) {
        return;
    }

    int left_idx = 0;
    int right_idx = 0;
    if (periodic_boundary != 0) {
        if (iface == 0 || iface == nx) {
            left_idx = nx - 1;
            right_idx = 0;
        } else {
            left_idx = iface - 1;
            right_idx = iface;
        }
    } else {
        if (iface == 0) {
            left_idx = 0;
            right_idx = 0;
        } else if (iface == nx) {
            left_idx = nx - 1;
            right_idx = nx - 1;
        } else {
            left_idx = iface - 1;
            right_idx = iface;
        }
    }

    fluxes[iface] = fp_cuda_hll_flux(cells[left_idx], cells[right_idx], gamma);
}

__global__ void fp_cuda_godunov_update_kernel(const fp::Conservative1D* cells,
                                              const fp::Conservative1D* fluxes,
                                              fp::Conservative1D* next_cells,
                                              int nx,
                                              double dt_over_dx) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nx) {
        return;
    }

    next_cells[i].rho = cells[i].rho - dt_over_dx * (fluxes[i + 1].rho - fluxes[i].rho);
    next_cells[i].mom = cells[i].mom - dt_over_dx * (fluxes[i + 1].mom - fluxes[i].mom);
    next_cells[i].energy = cells[i].energy - dt_over_dx * (fluxes[i + 1].energy - fluxes[i].energy);
}

__global__ void fp_cuda_fill_extended_kernel(const fp::Conservative1D* cells,
                                             fp::Conservative1D* extended,
                                             int nx,
                                             int periodic_boundary) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nx + 4) {
        return;
    }

    if (idx >= 2 && idx <= nx + 1) {
        extended[idx] = cells[idx - 2];
        return;
    }

    if (periodic_boundary != 0) {
        if (idx == 0) {
            extended[idx] = cells[nx - 2];
        } else if (idx == 1) {
            extended[idx] = cells[nx - 1];
        } else if (idx == nx + 2) {
            extended[idx] = cells[0];
        } else if (idx == nx + 3) {
            extended[idx] = cells[1];
        }
    } else {
        if (idx == 0 || idx == 1) {
            extended[idx] = cells[0];
        } else {
            extended[idx] = cells[nx - 1];
        }
    }
}

__global__ void fp_cuda_muscl_slopes_kernel(const fp::Conservative1D* extended,
                                            fp::Conservative1D* slopes,
                                            int nx) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if (i >= nx + 3) {
        return;
    }

    const fp::Conservative1D left_diff = {extended[i].rho - extended[i - 1].rho,
                                          extended[i].mom - extended[i - 1].mom,
                                          extended[i].energy - extended[i - 1].energy};
    const fp::Conservative1D right_diff = {extended[i + 1].rho - extended[i].rho,
                                           extended[i + 1].mom - extended[i].mom,
                                           extended[i + 1].energy - extended[i].energy};
    slopes[i] = fp_cuda_limited_slope(left_diff, right_diff);
}

__global__ void fp_cuda_muscl_predict_flux_kernel(const fp::Conservative1D* extended,
                                                  const fp::Conservative1D* slopes,
                                                  fp::Conservative1D* left_states,
                                                  fp::Conservative1D* right_states,
                                                  fp::Conservative1D* fluxes,
                                                  int nx,
                                                  double gamma,
                                                  double dt_over_dx) {
    const int iface = blockIdx.x * blockDim.x + threadIdx.x;
    if (iface > nx) {
        return;
    }

    const int cell_left = iface + 1;
    const int cell_right = iface + 2;

    const fp::Conservative1D u_left = extended[cell_left];
    const fp::Conservative1D u_right = extended[cell_right];
    const fp::Conservative1D slope_left = slopes[cell_left];
    const fp::Conservative1D slope_right = slopes[cell_right];

    const fp::Conservative1D u_left_plus = {u_left.rho + 0.5 * slope_left.rho,
                                            u_left.mom + 0.5 * slope_left.mom,
                                            u_left.energy + 0.5 * slope_left.energy};
    const fp::Conservative1D u_left_minus = {u_left.rho - 0.5 * slope_left.rho,
                                             u_left.mom - 0.5 * slope_left.mom,
                                             u_left.energy - 0.5 * slope_left.energy};
    const fp::Conservative1D u_right_plus = {u_right.rho + 0.5 * slope_right.rho,
                                             u_right.mom + 0.5 * slope_right.mom,
                                             u_right.energy + 0.5 * slope_right.energy};
    const fp::Conservative1D u_right_minus = {u_right.rho - 0.5 * slope_right.rho,
                                              u_right.mom - 0.5 * slope_right.mom,
                                              u_right.energy - 0.5 * slope_right.energy};

    const fp::Conservative1D left_flux_plus = fp_cuda_physical_flux(u_left_plus, gamma);
    const fp::Conservative1D left_flux_minus = fp_cuda_physical_flux(u_left_minus, gamma);
    const fp::Conservative1D right_flux_plus = fp_cuda_physical_flux(u_right_plus, gamma);
    const fp::Conservative1D right_flux_minus = fp_cuda_physical_flux(u_right_minus, gamma);

    left_states[iface] = {u_left_plus.rho - 0.5 * dt_over_dx * (left_flux_plus.rho - left_flux_minus.rho),
                          u_left_plus.mom - 0.5 * dt_over_dx * (left_flux_plus.mom - left_flux_minus.mom),
                          u_left_plus.energy -
                              0.5 * dt_over_dx * (left_flux_plus.energy - left_flux_minus.energy)};
    right_states[iface] = {u_right_minus.rho - 0.5 * dt_over_dx * (right_flux_plus.rho - right_flux_minus.rho),
                           u_right_minus.mom - 0.5 * dt_over_dx * (right_flux_plus.mom - right_flux_minus.mom),
                           u_right_minus.energy -
                               0.5 * dt_over_dx * (right_flux_plus.energy - right_flux_minus.energy)};
    fluxes[iface] = fp_cuda_hll_flux(left_states[iface], right_states[iface], gamma);
}

__global__ void fp_cuda_weno2_rhs_1d_kernel(const fp::Conservative1D* cells,
                                            fp::Conservative1D* rhs,
                                            int nx,
                                            double gamma,
                                            double inv_dx,
                                            int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nx) {
        return;
    }

    const fp::Conservative1D um2 = fp_cuda_sample_cell_1d(cells, nx, i - 2, periodic_boundary);
    const fp::Conservative1D um1 = fp_cuda_sample_cell_1d(cells, nx, i - 1, periodic_boundary);
    const fp::Conservative1D u0 = fp_cuda_sample_cell_1d(cells, nx, i, periodic_boundary);
    const fp::Conservative1D up1 = fp_cuda_sample_cell_1d(cells, nx, i + 1, periodic_boundary);
    const fp::Conservative1D up2 = fp_cuda_sample_cell_1d(cells, nx, i + 2, periodic_boundary);

    const fp::Conservative1D left_plus = fp_cuda_weno2_left_state_1d(um1, u0, up1);
    const fp::Conservative1D right_plus = fp_cuda_weno2_right_state_1d(u0, up1, up2);
    const fp::Conservative1D flux_plus = fp_cuda_hll_flux(left_plus, right_plus, gamma);

    const fp::Conservative1D left_minus = fp_cuda_weno2_left_state_1d(um2, um1, u0);
    const fp::Conservative1D right_minus = fp_cuda_weno2_right_state_1d(um1, u0, up1);
    const fp::Conservative1D flux_minus = fp_cuda_hll_flux(left_minus, right_minus, gamma);

    rhs[i].rho = -inv_dx * (flux_plus.rho - flux_minus.rho);
    rhs[i].mom = -inv_dx * (flux_plus.mom - flux_minus.mom);
    rhs[i].energy = -inv_dx * (flux_plus.energy - flux_minus.energy);
}

__global__ void fp_cuda_weno2_scalar_rhs_1d_kernel(const double* rho,
                                                   double* rhs,
                                                   int nx,
                                                   double advection_velocity,
                                                   double inv_dx,
                                                   int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nx) {
        return;
    }

    const double um2 = fp_cuda_sample_scalar_1d(rho, nx, i - 2, periodic_boundary);
    const double um1 = fp_cuda_sample_scalar_1d(rho, nx, i - 1, periodic_boundary);
    const double u0 = fp_cuda_sample_scalar_1d(rho, nx, i, periodic_boundary);
    const double up1 = fp_cuda_sample_scalar_1d(rho, nx, i + 1, periodic_boundary);
    const double up2 = fp_cuda_sample_scalar_1d(rho, nx, i + 2, periodic_boundary);

    const double left_plus = fp_cuda_weno2_left_component(um1, u0, up1);
    const double right_plus = fp_cuda_weno2_right_component(u0, up1, up2);
    const double flux_plus = (advection_velocity >= 0.0) ? advection_velocity * left_plus
                                                         : advection_velocity * right_plus;

    const double left_minus = fp_cuda_weno2_left_component(um2, um1, u0);
    const double right_minus = fp_cuda_weno2_right_component(um1, u0, up1);
    const double flux_minus = (advection_velocity >= 0.0) ? advection_velocity * left_minus
                                                          : advection_velocity * right_minus;

    rhs[i] = -inv_dx * (flux_plus - flux_minus);
}

__global__ void fp_cuda_add_scaled_1d_kernel(const fp::Conservative1D* base,
                                             const fp::Conservative1D* delta,
                                             fp::Conservative1D* out,
                                             int n,
                                             double scale) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    out[idx].rho = base[idx].rho + scale * delta[idx].rho;
    out[idx].mom = base[idx].mom + scale * delta[idx].mom;
    out[idx].energy = base[idx].energy + scale * delta[idx].energy;
}

__global__ void fp_cuda_rk2_finalize_1d_kernel(const fp::Conservative1D* field,
                                               const fp::Conservative1D* stage1,
                                               const fp::Conservative1D* rhs2,
                                               fp::Conservative1D* out,
                                               int n,
                                               double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    out[idx].rho = 0.5 * field[idx].rho + 0.5 * (stage1[idx].rho + dt * rhs2[idx].rho);
    out[idx].mom = 0.5 * field[idx].mom + 0.5 * (stage1[idx].mom + dt * rhs2[idx].mom);
    out[idx].energy = 0.5 * field[idx].energy + 0.5 * (stage1[idx].energy + dt * rhs2[idx].energy);
}

__global__ void fp_cuda_add_scaled_scalar_kernel(const double* base,
                                                 const double* delta,
                                                 double* out,
                                                 int n,
                                                 double scale) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    out[idx] = base[idx] + scale * delta[idx];
}

__global__ void fp_cuda_rk2_finalize_scalar_kernel(const double* field,
                                                   const double* stage1,
                                                   const double* rhs2,
                                                   double* out,
                                                   int n,
                                                   double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    out[idx] = 0.5 * field[idx] + 0.5 * (stage1[idx] + dt * rhs2[idx]);
}

__global__ void fp_cuda_godunov_flux_x_2d_kernel(const fp::Conservative2D* cells,
                                                 fp::Conservative2D* fx,
                                                 int nx,
                                                 int ny,
                                                 double gamma,
                                                 int periodic_boundary) {
    const int iface = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (iface > nx || j >= ny) {
        return;
    }

    int left_i = 0;
    int right_i = 0;
    if (periodic_boundary != 0) {
        if (iface == 0 || iface == nx) {
            left_i = nx - 1;
            right_i = 0;
        } else {
            left_i = iface - 1;
            right_i = iface;
        }
    } else {
        if (iface == 0) {
            left_i = 0;
            right_i = 0;
        } else if (iface == nx) {
            left_i = nx - 1;
            right_i = nx - 1;
        } else {
            left_i = iface - 1;
            right_i = iface;
        }
    }

    const std::size_t left_idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + left_i;
    const std::size_t right_idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + right_i;
    const std::size_t flux_idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx + 1) + iface;
    fx[flux_idx] = fp_cuda_rusanov_flux_x_2d(cells[left_idx], cells[right_idx], gamma);
}

__global__ void fp_cuda_godunov_flux_y_2d_kernel(const fp::Conservative2D* cells,
                                                 fp::Conservative2D* fy,
                                                 int nx,
                                                 int ny,
                                                 double gamma,
                                                 int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int iface = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || iface > ny) {
        return;
    }

    int bottom_j = 0;
    int top_j = 0;
    if (periodic_boundary != 0) {
        if (iface == 0 || iface == ny) {
            bottom_j = ny - 1;
            top_j = 0;
        } else {
            bottom_j = iface - 1;
            top_j = iface;
        }
    } else {
        if (iface == 0) {
            bottom_j = 0;
            top_j = 0;
        } else if (iface == ny) {
            bottom_j = ny - 1;
            top_j = ny - 1;
        } else {
            bottom_j = iface - 1;
            top_j = iface;
        }
    }

    const std::size_t bottom_idx = static_cast<std::size_t>(bottom_j) * static_cast<std::size_t>(nx) + i;
    const std::size_t top_idx = static_cast<std::size_t>(top_j) * static_cast<std::size_t>(nx) + i;
    const std::size_t flux_idx = static_cast<std::size_t>(iface) * static_cast<std::size_t>(nx) + i;
    fy[flux_idx] = fp_cuda_rusanov_flux_y_2d(cells[bottom_idx], cells[top_idx], gamma);
}

__global__ void fp_cuda_godunov_update_2d_kernel(const fp::Conservative2D* cells,
                                                 const fp::Conservative2D* fx,
                                                 const fp::Conservative2D* fy,
                                                 fp::Conservative2D* next_cells,
                                                 int nx,
                                                 int ny,
                                                 double dt_over_dx,
                                                 double dt_over_dy) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) {
        return;
    }

    const std::size_t idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + i;
    const std::size_t fx_l = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx + 1) + i;
    const std::size_t fx_r = fx_l + 1;
    const std::size_t fy_b = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + i;
    const std::size_t fy_t = static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(nx) + i;

    next_cells[idx].rho =
        cells[idx].rho - dt_over_dx * (fx[fx_r].rho - fx[fx_l].rho) - dt_over_dy * (fy[fy_t].rho - fy[fy_b].rho);
    next_cells[idx].mom_x = cells[idx].mom_x - dt_over_dx * (fx[fx_r].mom_x - fx[fx_l].mom_x) -
                            dt_over_dy * (fy[fy_t].mom_x - fy[fy_b].mom_x);
    next_cells[idx].mom_y = cells[idx].mom_y - dt_over_dx * (fx[fx_r].mom_y - fx[fx_l].mom_y) -
                            dt_over_dy * (fy[fy_t].mom_y - fy[fy_b].mom_y);
    next_cells[idx].energy = cells[idx].energy - dt_over_dx * (fx[fx_r].energy - fx[fx_l].energy) -
                             dt_over_dy * (fy[fy_t].energy - fy[fy_b].energy);
}

__global__ void fp_cuda_muscl_rhs_2d_kernel(const fp::Conservative2D* cells,
                                            fp::Conservative2D* rhs,
                                            int nx,
                                            int ny,
                                            double gamma,
                                            double inv_dx,
                                            double inv_dy,
                                            int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) {
        return;
    }

    const fp::Conservative2D left_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    const fp::Conservative2D center_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D right_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    const fp::Conservative2D right_right_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i + 2, j, periodic_boundary);

    const fp::Conservative2D slope_center_x = fp_cuda_limited_slope_x_2d(left_cell, center_cell, right_cell);
    const fp::Conservative2D slope_right_x = fp_cuda_limited_slope_x_2d(center_cell, right_cell, right_right_cell);

    const fp::Conservative2D left_state_x = {center_cell.rho + 0.5 * slope_center_x.rho,
                                             center_cell.mom_x + 0.5 * slope_center_x.mom_x,
                                             center_cell.mom_y + 0.5 * slope_center_x.mom_y,
                                             center_cell.energy + 0.5 * slope_center_x.energy};
    const fp::Conservative2D right_state_x = {right_cell.rho - 0.5 * slope_right_x.rho,
                                              right_cell.mom_x - 0.5 * slope_right_x.mom_x,
                                              right_cell.mom_y - 0.5 * slope_right_x.mom_y,
                                              right_cell.energy - 0.5 * slope_right_x.energy};
    const fp::Conservative2D flux_x_plus = fp_cuda_rusanov_flux_x_2d(left_state_x, right_state_x, gamma);

    const fp::Conservative2D left_left_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i - 2, j, periodic_boundary);
    const fp::Conservative2D slope_left_x = fp_cuda_limited_slope_x_2d(left_left_cell, left_cell, center_cell);
    const fp::Conservative2D left_state_xm = {left_cell.rho + 0.5 * slope_left_x.rho,
                                              left_cell.mom_x + 0.5 * slope_left_x.mom_x,
                                              left_cell.mom_y + 0.5 * slope_left_x.mom_y,
                                              left_cell.energy + 0.5 * slope_left_x.energy};
    const fp::Conservative2D right_state_xm = {center_cell.rho - 0.5 * slope_center_x.rho,
                                               center_cell.mom_x - 0.5 * slope_center_x.mom_x,
                                               center_cell.mom_y - 0.5 * slope_center_x.mom_y,
                                               center_cell.energy - 0.5 * slope_center_x.energy};
    const fp::Conservative2D flux_x_minus = fp_cuda_rusanov_flux_x_2d(left_state_xm, right_state_xm, gamma);

    const fp::Conservative2D bottom_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    const fp::Conservative2D top_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    const fp::Conservative2D top_top_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 2, periodic_boundary);

    const fp::Conservative2D slope_center_y = fp_cuda_limited_slope_y_2d(bottom_cell, center_cell, top_cell);
    const fp::Conservative2D slope_top_y = fp_cuda_limited_slope_y_2d(center_cell, top_cell, top_top_cell);

    const fp::Conservative2D bottom_state_y = {center_cell.rho + 0.5 * slope_center_y.rho,
                                               center_cell.mom_x + 0.5 * slope_center_y.mom_x,
                                               center_cell.mom_y + 0.5 * slope_center_y.mom_y,
                                               center_cell.energy + 0.5 * slope_center_y.energy};
    const fp::Conservative2D top_state_y = {top_cell.rho - 0.5 * slope_top_y.rho,
                                            top_cell.mom_x - 0.5 * slope_top_y.mom_x,
                                            top_cell.mom_y - 0.5 * slope_top_y.mom_y,
                                            top_cell.energy - 0.5 * slope_top_y.energy};
    const fp::Conservative2D flux_y_plus = fp_cuda_rusanov_flux_y_2d(bottom_state_y, top_state_y, gamma);

    const fp::Conservative2D bottom_bottom_cell = fp_cuda_sample_cell_2d(cells, nx, ny, i, j - 2, periodic_boundary);
    const fp::Conservative2D slope_bottom_y = fp_cuda_limited_slope_y_2d(bottom_bottom_cell, bottom_cell, center_cell);
    const fp::Conservative2D bottom_state_ym = {bottom_cell.rho + 0.5 * slope_bottom_y.rho,
                                                bottom_cell.mom_x + 0.5 * slope_bottom_y.mom_x,
                                                bottom_cell.mom_y + 0.5 * slope_bottom_y.mom_y,
                                                bottom_cell.energy + 0.5 * slope_bottom_y.energy};
    const fp::Conservative2D top_state_ym = {center_cell.rho - 0.5 * slope_center_y.rho,
                                             center_cell.mom_x - 0.5 * slope_center_y.mom_x,
                                             center_cell.mom_y - 0.5 * slope_center_y.mom_y,
                                             center_cell.energy - 0.5 * slope_center_y.energy};
    const fp::Conservative2D flux_y_minus = fp_cuda_rusanov_flux_y_2d(bottom_state_ym, top_state_ym, gamma);

    const std::size_t idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i);
    rhs[idx].rho = -inv_dx * (flux_x_plus.rho - flux_x_minus.rho) - inv_dy * (flux_y_plus.rho - flux_y_minus.rho);
    rhs[idx].mom_x =
        -inv_dx * (flux_x_plus.mom_x - flux_x_minus.mom_x) - inv_dy * (flux_y_plus.mom_x - flux_y_minus.mom_x);
    rhs[idx].mom_y =
        -inv_dx * (flux_x_plus.mom_y - flux_x_minus.mom_y) - inv_dy * (flux_y_plus.mom_y - flux_y_minus.mom_y);
    rhs[idx].energy = -inv_dx * (flux_x_plus.energy - flux_x_minus.energy) -
                      inv_dy * (flux_y_plus.energy - flux_y_minus.energy);
}

__global__ void fp_cuda_weno2_rhs_2d_kernel(const fp::Conservative2D* cells,
                                            fp::Conservative2D* rhs,
                                            int nx,
                                            int ny,
                                            double gamma,
                                            double inv_dx,
                                            double inv_dy,
                                            int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) {
        return;
    }

    const fp::Conservative2D left_plus =
        fp_cuda_weno2_left_state_x_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D right_plus =
        fp_cuda_weno2_right_state_x_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D flux_x_plus = fp_cuda_rusanov_flux_x_2d(left_plus, right_plus, gamma);

    const fp::Conservative2D left_minus =
        fp_cuda_weno2_left_state_x_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    const fp::Conservative2D right_minus =
        fp_cuda_weno2_right_state_x_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    const fp::Conservative2D flux_x_minus = fp_cuda_rusanov_flux_x_2d(left_minus, right_minus, gamma);

    const fp::Conservative2D bottom_plus =
        fp_cuda_weno2_left_state_y_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D top_plus =
        fp_cuda_weno2_right_state_y_2d(cells, nx, ny, i, j, periodic_boundary);
    const fp::Conservative2D flux_y_plus = fp_cuda_rusanov_flux_y_2d(bottom_plus, top_plus, gamma);

    const fp::Conservative2D bottom_minus =
        fp_cuda_weno2_left_state_y_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    const fp::Conservative2D top_minus =
        fp_cuda_weno2_right_state_y_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    const fp::Conservative2D flux_y_minus = fp_cuda_rusanov_flux_y_2d(bottom_minus, top_minus, gamma);

    const std::size_t idx = static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i);
    rhs[idx].rho = -inv_dx * (flux_x_plus.rho - flux_x_minus.rho) - inv_dy * (flux_y_plus.rho - flux_y_minus.rho);
    rhs[idx].mom_x =
        -inv_dx * (flux_x_plus.mom_x - flux_x_minus.mom_x) - inv_dy * (flux_y_plus.mom_x - flux_y_minus.mom_x);
    rhs[idx].mom_y =
        -inv_dx * (flux_x_plus.mom_y - flux_x_minus.mom_y) - inv_dy * (flux_y_plus.mom_y - flux_y_minus.mom_y);
    rhs[idx].energy = -inv_dx * (flux_x_plus.energy - flux_x_minus.energy) -
                      inv_dy * (flux_y_plus.energy - flux_y_minus.energy);
}

__global__ void fp_cuda_weno3_rhs_2d_kernel(const fp::Conservative2D* cells,
                                            fp::Conservative2D* rhs,
                                            int nx,
                                            int ny,
                                            double gamma,
                                            double inv_dx,
                                            double inv_dy,
                                            int periodic_boundary) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) {
        return;
    }

    fp::Conservative2D left_plus =
        fp_cuda_weno3_left_state_x_2d(cells, nx, ny, i, j, periodic_boundary);
    fp::Conservative2D right_plus =
        fp_cuda_weno3_right_state_x_2d(cells, nx, ny, i, j, periodic_boundary);
    if (!fp_cuda_is_admissible_2d(left_plus, gamma)) {
        left_plus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    }
    if (!fp_cuda_is_admissible_2d(right_plus, gamma)) {
        right_plus = fp_cuda_sample_cell_2d(cells, nx, ny, i + 1, j, periodic_boundary);
    }
    const fp::Conservative2D flux_x_plus =
        fp_cuda_rusanov_flux_x_2d(left_plus, right_plus, gamma);

    fp::Conservative2D left_minus =
        fp_cuda_weno3_left_state_x_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    fp::Conservative2D right_minus =
        fp_cuda_weno3_right_state_x_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    if (!fp_cuda_is_admissible_2d(left_minus, gamma)) {
        left_minus = fp_cuda_sample_cell_2d(cells, nx, ny, i - 1, j, periodic_boundary);
    }
    if (!fp_cuda_is_admissible_2d(right_minus, gamma)) {
        right_minus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    }
    const fp::Conservative2D flux_x_minus =
        fp_cuda_rusanov_flux_x_2d(left_minus, right_minus, gamma);

    fp::Conservative2D bottom_plus =
        fp_cuda_weno3_left_state_y_2d(cells, nx, ny, i, j, periodic_boundary);
    fp::Conservative2D top_plus =
        fp_cuda_weno3_right_state_y_2d(cells, nx, ny, i, j, periodic_boundary);
    if (!fp_cuda_is_admissible_2d(bottom_plus, gamma)) {
        bottom_plus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    }
    if (!fp_cuda_is_admissible_2d(top_plus, gamma)) {
        top_plus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j + 1, periodic_boundary);
    }
    const fp::Conservative2D flux_y_plus =
        fp_cuda_rusanov_flux_y_2d(bottom_plus, top_plus, gamma);

    fp::Conservative2D bottom_minus =
        fp_cuda_weno3_left_state_y_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    fp::Conservative2D top_minus =
        fp_cuda_weno3_right_state_y_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    if (!fp_cuda_is_admissible_2d(bottom_minus, gamma)) {
        bottom_minus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j - 1, periodic_boundary);
    }
    if (!fp_cuda_is_admissible_2d(top_minus, gamma)) {
        top_minus = fp_cuda_sample_cell_2d(cells, nx, ny, i, j, periodic_boundary);
    }
    const fp::Conservative2D flux_y_minus =
        fp_cuda_rusanov_flux_y_2d(bottom_minus, top_minus, gamma);

    const std::size_t idx =
        static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) + static_cast<std::size_t>(i);
    rhs[idx].rho =
        -inv_dx * (flux_x_plus.rho - flux_x_minus.rho) -
        inv_dy * (flux_y_plus.rho - flux_y_minus.rho);
    rhs[idx].mom_x =
        -inv_dx * (flux_x_plus.mom_x - flux_x_minus.mom_x) -
        inv_dy * (flux_y_plus.mom_x - flux_y_minus.mom_x);
    rhs[idx].mom_y =
        -inv_dx * (flux_x_plus.mom_y - flux_x_minus.mom_y) -
        inv_dy * (flux_y_plus.mom_y - flux_y_minus.mom_y);
    rhs[idx].energy =
        -inv_dx * (flux_x_plus.energy - flux_x_minus.energy) -
        inv_dy * (flux_y_plus.energy - flux_y_minus.energy);
}

__device__ double fp_cuda_atomic_max_nonnegative_double(double* address, double value) {
    auto* address_as_ull = reinterpret_cast<unsigned long long int*>(address);
    unsigned long long int old = *address_as_ull;
    while (__longlong_as_double(old) < value) {
        const unsigned long long int assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(value));
        if (old == assumed) {
            break;
        }
    }
    return __longlong_as_double(old);
}

__global__ void fp_cuda_max_signal_speed_2d_kernel(const fp::Conservative2D* cells,
                                                   int n,
                                                   double gamma,
                                                   double* max_x,
                                                   double* max_y) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    fp_cuda_atomic_max_nonnegative_double(
        max_x, fp_cuda_max_signal_speed_x_2d(cells[idx], gamma));
    fp_cuda_atomic_max_nonnegative_double(
        max_y, fp_cuda_max_signal_speed_y_2d(cells[idx], gamma));
}

__global__ void fp_cuda_add_scaled_2d_kernel(const fp::Conservative2D* base,
                                             const fp::Conservative2D* delta,
                                             fp::Conservative2D* out,
                                             int n,
                                             double scale) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    out[idx].rho = base[idx].rho + scale * delta[idx].rho;
    out[idx].mom_x = base[idx].mom_x + scale * delta[idx].mom_x;
    out[idx].mom_y = base[idx].mom_y + scale * delta[idx].mom_y;
    out[idx].energy = base[idx].energy + scale * delta[idx].energy;
}

__global__ void fp_cuda_rk3_stage2_2d_kernel(const fp::Conservative2D* field,
                                             const fp::Conservative2D* stage1,
                                             const fp::Conservative2D* rhs2,
                                             fp::Conservative2D* stage2,
                                             int n,
                                             double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    stage2[idx].rho =
        0.75 * field[idx].rho + 0.25 * (stage1[idx].rho + dt * rhs2[idx].rho);
    stage2[idx].mom_x =
        0.75 * field[idx].mom_x + 0.25 * (stage1[idx].mom_x + dt * rhs2[idx].mom_x);
    stage2[idx].mom_y =
        0.75 * field[idx].mom_y + 0.25 * (stage1[idx].mom_y + dt * rhs2[idx].mom_y);
    stage2[idx].energy =
        0.75 * field[idx].energy + 0.25 * (stage1[idx].energy + dt * rhs2[idx].energy);
}

__global__ void fp_cuda_rk3_finalize_2d_kernel(const fp::Conservative2D* field,
                                               const fp::Conservative2D* stage2,
                                               const fp::Conservative2D* rhs3,
                                               fp::Conservative2D* out,
                                               int n,
                                               double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    out[idx].rho =
        (1.0 / 3.0) * field[idx].rho +
        (2.0 / 3.0) * (stage2[idx].rho + dt * rhs3[idx].rho);
    out[idx].mom_x =
        (1.0 / 3.0) * field[idx].mom_x +
        (2.0 / 3.0) * (stage2[idx].mom_x + dt * rhs3[idx].mom_x);
    out[idx].mom_y =
        (1.0 / 3.0) * field[idx].mom_y +
        (2.0 / 3.0) * (stage2[idx].mom_y + dt * rhs3[idx].mom_y);
    out[idx].energy =
        (1.0 / 3.0) * field[idx].energy +
        (2.0 / 3.0) * (stage2[idx].energy + dt * rhs3[idx].energy);
}

__global__ void fp_cuda_rk2_finalize_2d_kernel(const fp::Conservative2D* field,
                                               const fp::Conservative2D* stage1,
                                               const fp::Conservative2D* rhs2,
                                               fp::Conservative2D* out,
                                               int n,
                                               double dt) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }

    out[idx].rho = 0.5 * field[idx].rho + 0.5 * (stage1[idx].rho + dt * rhs2[idx].rho);
    out[idx].mom_x = 0.5 * field[idx].mom_x + 0.5 * (stage1[idx].mom_x + dt * rhs2[idx].mom_x);
    out[idx].mom_y = 0.5 * field[idx].mom_y + 0.5 * (stage1[idx].mom_y + dt * rhs2[idx].mom_y);
    out[idx].energy = 0.5 * field[idx].energy + 0.5 * (stage1[idx].energy + dt * rhs2[idx].energy);
}

double fp_host_compute_time_step(const std::vector<fp::Conservative1D>& cells,
                                 const fp::RuntimeConfig& config,
                                 double dx) {
    double max_speed = 1.0e-12;

    for (int iface = 0; iface <= config.nx; ++iface) {
        int left_idx = 0;
        int right_idx = 0;
        if (config.boundary == "periodic") {
            if (iface == 0 || iface == config.nx) {
                left_idx = config.nx - 1;
                right_idx = 0;
            } else {
                left_idx = iface - 1;
                right_idx = iface;
            }
        } else {
            if (iface == 0) {
                left_idx = 0;
                right_idx = 0;
            } else if (iface == config.nx) {
                left_idx = config.nx - 1;
                right_idx = config.nx - 1;
            } else {
                left_idx = iface - 1;
                right_idx = iface;
            }
        }

        max_speed = std::max(max_speed, fp::max_wave_speed(cells[left_idx], cells[right_idx], config.gamma));
    }

    return config.cfl * dx / max_speed;
}

double fp_host_compute_time_step_scalar(const fp::RuntimeConfig& config, double dx) {
    const double speed = std::abs(config.advection_velocity);
    if (speed < 1.0e-14) {
        return config.final_time;
    }
    return config.cfl * dx / speed;
}

double fp_host_compute_time_step_2d(const fp::Field2D& field,
                                    const fp::Grid2D& grid,
                                    const fp::RuntimeConfig& config) {
    double max_speed = 1.0e-12;
    for (int j = 0; j < field.ny; ++j) {
        for (int i = 0; i < field.nx; ++i) {
            max_speed = std::max(max_speed, fp::max_signal_speed_x(field(i, j), config.gamma));
            max_speed = std::max(max_speed, fp::max_signal_speed_y(field(i, j), config.gamma));
        }
    }
    return config.cfl * std::min(fp::dx(grid), fp::dy(grid)) / max_speed;
}

}  // namespace

namespace fp {

bool cuda_backend_compiled() {
    return true;
}

bool cuda_backend_probe() {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0) {
        return false;
    }

    int* device_flag = nullptr;
    int host_flag = 0;
    if (cudaMalloc(&device_flag, sizeof(int)) != cudaSuccess) {
        return false;
    }

    if (cudaMemset(device_flag, 0, sizeof(int)) != cudaSuccess) {
        cudaFree(device_flag);
        return false;
    }

    fp_cuda_smoke_kernel<<<1, 1>>>(device_flag);
    if (cudaDeviceSynchronize() != cudaSuccess) {
        cudaFree(device_flag);
        return false;
    }

    if (cudaMemcpy(&host_flag, device_flag, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) {
        cudaFree(device_flag);
        return false;
    }

    cudaFree(device_flag);
    return host_flag == 1;
}

std::string cuda_backend_summary() {
    int device_count = 0;
    cudaError_t status = cudaGetDeviceCount(&device_count);
    if (status != cudaSuccess) {
        return "CUDA compiled, but device query failed.";
    }

    std::ostringstream oss;
    oss << "CUDA compiled. Detected device count: " << device_count;
    return oss.str();
}

SimulationResult1D cuda_copy_1d_state(const RuntimeConfig& config,
                                      const std::vector<Conservative1D>& initial_state) {
    if (static_cast<int>(initial_state.size()) != config.nx) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    Conservative1D* d_in = nullptr;
    Conservative1D* d_out = nullptr;
    const std::size_t bytes = sizeof(Conservative1D) * initial_state.size();

    if (cudaMalloc(&d_in, bytes) != cudaSuccess || cudaMalloc(&d_out, bytes) != cudaSuccess) {
        if (d_in) cudaFree(d_in);
        if (d_out) cudaFree(d_out);
        throw std::runtime_error("cudaMalloc failed in cuda_copy_1d_state.");
    }

    if (cudaMemcpy(d_in, initial_state.data(), bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_in);
        cudaFree(d_out);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_copy_1d_state.");
    }

    const int block_size = 256;
    const int grid_size = (config.nx + block_size - 1) / block_size;
    fp_cuda_copy_state_kernel<<<grid_size, block_size>>>(d_in, d_out, config.nx);
    if (cudaDeviceSynchronize() != cudaSuccess) {
        cudaFree(d_in);
        cudaFree(d_out);
        throw std::runtime_error("cudaDeviceSynchronize failed in cuda_copy_1d_state.");
    }

    SimulationResult1D result;
    result.cells.resize(initial_state.size());
    if (cudaMemcpy(result.cells.data(), d_out, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        cudaFree(d_in);
        cudaFree(d_out);
        throw std::runtime_error("cudaMemcpy D2H failed in cuda_copy_1d_state.");
    }

    cudaFree(d_in);
    cudaFree(d_out);

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    result.x.resize(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }
    result.final_time = 0.0;
    result.steps = 0;
    return result;
}

SimulationResult1D cuda_run_godunov_1d(const RuntimeConfig& config,
                                       const std::vector<Conservative1D>& initial_state) {
    if (static_cast<int>(initial_state.size()) != config.nx) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const std::size_t cell_bytes = sizeof(Conservative1D) * initial_state.size();
    const std::size_t flux_bytes = sizeof(Conservative1D) * static_cast<std::size_t>(config.nx + 1);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const int block_size = 256;
    const int cell_grid_size = (config.nx + block_size - 1) / block_size;
    const int flux_grid_size = (config.nx + 1 + block_size - 1) / block_size;
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative1D* d_cells = nullptr;
    Conservative1D* d_next = nullptr;
    Conservative1D* d_fluxes = nullptr;

    if (cudaMalloc(&d_cells, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_fluxes, flux_bytes) != cudaSuccess) {
        if (d_cells) cudaFree(d_cells);
        if (d_next) cudaFree(d_next);
        if (d_fluxes) cudaFree(d_fluxes);
        throw std::runtime_error("cudaMalloc failed in cuda_run_godunov_1d.");
    }

    if (cudaMemcpy(d_cells, initial_state.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_cells);
        cudaFree(d_next);
        cudaFree(d_fluxes);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_godunov_1d.");
    }

    std::vector<Conservative1D> host_cells = initial_state;
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemcpy(host_cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            cudaFree(d_cells);
            cudaFree(d_next);
            cudaFree(d_fluxes);
            throw std::runtime_error("cudaMemcpy D2H state failed in cuda_run_godunov_1d.");
        }

        double dt = fp_host_compute_time_step(host_cells, config, dx);
        dt = std::min(dt, config.final_time - time);

        fp_cuda_godunov_flux_kernel<<<flux_grid_size, block_size>>>(d_cells, d_fluxes, config.nx, config.gamma,
                                                                    periodic_boundary);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            cudaFree(d_cells);
            cudaFree(d_next);
            cudaFree(d_fluxes);
            throw std::runtime_error("Flux kernel failed in cuda_run_godunov_1d.");
        }

        fp_cuda_godunov_update_kernel<<<cell_grid_size, block_size>>>(d_cells, d_fluxes, d_next, config.nx, dt / dx);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            cudaFree(d_cells);
            cudaFree(d_next);
            cudaFree(d_fluxes);
            throw std::runtime_error("Update kernel failed in cuda_run_godunov_1d.");
        }

        std::swap(d_cells, d_next);
        time += dt;
        ++step;
    }

    SimulationResult1D result;
    result.cells.resize(initial_state.size());
    if (cudaMemcpy(result.cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        cudaFree(d_cells);
        cudaFree(d_next);
        cudaFree(d_fluxes);
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_godunov_1d.");
    }

    cudaFree(d_cells);
    cudaFree(d_next);
    cudaFree(d_fluxes);

    result.x.resize(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }
    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult1D cuda_run_muscl_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state) {
    if (static_cast<int>(initial_state.size()) != config.nx) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const std::size_t cell_bytes = sizeof(Conservative1D) * initial_state.size();
    const std::size_t extended_bytes = sizeof(Conservative1D) * static_cast<std::size_t>(config.nx + 4);
    const std::size_t slope_bytes = sizeof(Conservative1D) * static_cast<std::size_t>(config.nx + 4);
    const std::size_t iface_bytes = sizeof(Conservative1D) * static_cast<std::size_t>(config.nx + 1);
    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const int block_size = 256;
    const int cell_grid_size = (config.nx + block_size - 1) / block_size;
    const int extended_grid_size = (config.nx + 4 + block_size - 1) / block_size;
    const int slope_grid_size = (config.nx + 2 + block_size - 1) / block_size;
    const int flux_grid_size = (config.nx + 1 + block_size - 1) / block_size;
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative1D* d_cells = nullptr;
    Conservative1D* d_next = nullptr;
    Conservative1D* d_extended = nullptr;
    Conservative1D* d_slopes = nullptr;
    Conservative1D* d_left_states = nullptr;
    Conservative1D* d_right_states = nullptr;
    Conservative1D* d_fluxes = nullptr;

    if (cudaMalloc(&d_cells, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_extended, extended_bytes) != cudaSuccess ||
        cudaMalloc(&d_slopes, slope_bytes) != cudaSuccess ||
        cudaMalloc(&d_left_states, iface_bytes) != cudaSuccess ||
        cudaMalloc(&d_right_states, iface_bytes) != cudaSuccess ||
        cudaMalloc(&d_fluxes, iface_bytes) != cudaSuccess) {
        if (d_cells) cudaFree(d_cells);
        if (d_next) cudaFree(d_next);
        if (d_extended) cudaFree(d_extended);
        if (d_slopes) cudaFree(d_slopes);
        if (d_left_states) cudaFree(d_left_states);
        if (d_right_states) cudaFree(d_right_states);
        if (d_fluxes) cudaFree(d_fluxes);
        throw std::runtime_error("cudaMalloc failed in cuda_run_muscl_1d.");
    }

    if (cudaMemcpy(d_cells, initial_state.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_cells);
        cudaFree(d_next);
        cudaFree(d_extended);
        cudaFree(d_slopes);
        cudaFree(d_left_states);
        cudaFree(d_right_states);
        cudaFree(d_fluxes);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_muscl_1d.");
    }

    std::vector<Conservative1D> host_cells = initial_state;
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemcpy(host_cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            throw std::runtime_error("cudaMemcpy D2H state failed in cuda_run_muscl_1d.");
        }

        double dt = fp_host_compute_time_step(host_cells, config, dx);
        dt = std::min(dt, config.final_time - time);

        fp_cuda_fill_extended_kernel<<<extended_grid_size, block_size>>>(d_cells, d_extended, config.nx,
                                                                         periodic_boundary);
        fp_cuda_muscl_slopes_kernel<<<slope_grid_size, block_size>>>(d_extended, d_slopes, config.nx);
        fp_cuda_muscl_predict_flux_kernel<<<flux_grid_size, block_size>>>(
            d_extended, d_slopes, d_left_states, d_right_states, d_fluxes, config.nx, config.gamma, dt / dx);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Reconstruction/flux kernels failed in cuda_run_muscl_1d.");
        }

        fp_cuda_godunov_update_kernel<<<cell_grid_size, block_size>>>(d_cells, d_fluxes, d_next, config.nx, dt / dx);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Update kernel failed in cuda_run_muscl_1d.");
        }

        std::swap(d_cells, d_next);
        time += dt;
        ++step;
    }

    SimulationResult1D result;
    result.cells.resize(initial_state.size());
    if (cudaMemcpy(result.cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_muscl_1d.");
    }

    cudaFree(d_cells);
    cudaFree(d_next);
    cudaFree(d_extended);
    cudaFree(d_slopes);
    cudaFree(d_left_states);
    cudaFree(d_right_states);
    cudaFree(d_fluxes);

    result.x.resize(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }
    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult1D cuda_run_weno2_1d(const RuntimeConfig& config,
                                     const std::vector<Conservative1D>& initial_state) {
    if (static_cast<int>(initial_state.size()) != config.nx) {
        throw std::runtime_error("Initial state size does not match config.nx.");
    }

    const double dx = (config.x_max - config.x_min) / static_cast<double>(config.nx);
    const int block_size = 256;
    const int grid_size = (config.nx + block_size - 1) / block_size;
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;
    const bool use_scalar_advection_path =
        (config.problem == "smooth_advection" || config.problem == "gaussian_advection") && periodic_boundary != 0;

    SimulationResult1D result;
    result.cells.resize(initial_state.size());
    result.x.resize(config.nx);
    for (int i = 0; i < config.nx; ++i) {
        result.x[i] = config.x_min + (static_cast<double>(i) + 0.5) * dx;
    }

    if (use_scalar_advection_path) {
        const std::size_t scalar_bytes = sizeof(double) * static_cast<std::size_t>(config.nx);
        std::vector<double> host_rho(static_cast<std::size_t>(config.nx), 0.0);
        for (int i = 0; i < config.nx; ++i) {
            host_rho[static_cast<std::size_t>(i)] = initial_state[static_cast<std::size_t>(i)].rho;
        }

        double* d_rho = nullptr;
        double* d_rhs1 = nullptr;
        double* d_stage1 = nullptr;
        double* d_rhs2 = nullptr;
        double* d_next = nullptr;

        if (cudaMalloc(&d_rho, scalar_bytes) != cudaSuccess ||
            cudaMalloc(&d_rhs1, scalar_bytes) != cudaSuccess ||
            cudaMalloc(&d_stage1, scalar_bytes) != cudaSuccess ||
            cudaMalloc(&d_rhs2, scalar_bytes) != cudaSuccess ||
            cudaMalloc(&d_next, scalar_bytes) != cudaSuccess) {
            if (d_rho) cudaFree(d_rho);
            if (d_rhs1) cudaFree(d_rhs1);
            if (d_stage1) cudaFree(d_stage1);
            if (d_rhs2) cudaFree(d_rhs2);
            if (d_next) cudaFree(d_next);
            throw std::runtime_error("cudaMalloc failed in cuda_run_weno2_1d scalar path.");
        }

        if (cudaMemcpy(d_rho, host_rho.data(), scalar_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
            cudaFree(d_rho);
            cudaFree(d_rhs1);
            cudaFree(d_stage1);
            cudaFree(d_rhs2);
            cudaFree(d_next);
            throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_weno2_1d scalar path.");
        }

        double time = 0.0;
        int step = 0;
        while (time < config.final_time && step < config.max_steps) {
            double dt = fp_host_compute_time_step_scalar(config, dx);
            dt = std::min(dt, config.final_time - time);

            fp_cuda_weno2_scalar_rhs_1d_kernel<<<grid_size, block_size>>>(d_rho, d_rhs1, config.nx,
                                                                           config.advection_velocity, 1.0 / dx,
                                                                           periodic_boundary);
            fp_cuda_add_scaled_scalar_kernel<<<grid_size, block_size>>>(d_rho, d_rhs1, d_stage1, config.nx, dt);
            if (cudaDeviceSynchronize() != cudaSuccess) {
                throw std::runtime_error("Stage 1 kernels failed in cuda_run_weno2_1d scalar path.");
            }

            fp_cuda_weno2_scalar_rhs_1d_kernel<<<grid_size, block_size>>>(d_stage1, d_rhs2, config.nx,
                                                                           config.advection_velocity, 1.0 / dx,
                                                                           periodic_boundary);
            fp_cuda_rk2_finalize_scalar_kernel<<<grid_size, block_size>>>(d_rho, d_stage1, d_rhs2, d_next, config.nx,
                                                                          dt);
            if (cudaDeviceSynchronize() != cudaSuccess) {
                throw std::runtime_error("Stage 2 kernels failed in cuda_run_weno2_1d scalar path.");
            }

            std::swap(d_rho, d_next);
            time += dt;
            ++step;
        }

        if (cudaMemcpy(host_rho.data(), d_rho, scalar_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_weno2_1d scalar path.");
        }

        cudaFree(d_rho);
        cudaFree(d_rhs1);
        cudaFree(d_stage1);
        cudaFree(d_rhs2);
        cudaFree(d_next);

        for (int i = 0; i < config.nx; ++i) {
            const Primitive1D state{host_rho[static_cast<std::size_t>(i)], config.advection_velocity,
                                    kBackgroundPressure};
            result.cells[static_cast<std::size_t>(i)] = primitive_to_conservative(state, config.gamma);
        }
        result.final_time = time;
        result.steps = step;
        return result;
    }

    const std::size_t cell_bytes = sizeof(Conservative1D) * initial_state.size();
    Conservative1D* d_cells = nullptr;
    Conservative1D* d_rhs1 = nullptr;
    Conservative1D* d_stage1 = nullptr;
    Conservative1D* d_rhs2 = nullptr;
    Conservative1D* d_next = nullptr;

    if (cudaMalloc(&d_cells, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_stage1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs2, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess) {
        if (d_cells) cudaFree(d_cells);
        if (d_rhs1) cudaFree(d_rhs1);
        if (d_stage1) cudaFree(d_stage1);
        if (d_rhs2) cudaFree(d_rhs2);
        if (d_next) cudaFree(d_next);
        throw std::runtime_error("cudaMalloc failed in cuda_run_weno2_1d.");
    }

    if (cudaMemcpy(d_cells, initial_state.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_cells);
        cudaFree(d_rhs1);
        cudaFree(d_stage1);
        cudaFree(d_rhs2);
        cudaFree(d_next);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_weno2_1d.");
    }

    std::vector<Conservative1D> host_cells = initial_state;
    double time = 0.0;
    int step = 0;

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemcpy(host_cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            throw std::runtime_error("cudaMemcpy D2H state failed in cuda_run_weno2_1d.");
        }

        double dt = fp_host_compute_time_step(host_cells, config, dx);
        dt = std::min(dt, config.final_time - time);

        fp_cuda_weno2_rhs_1d_kernel<<<grid_size, block_size>>>(d_cells, d_rhs1, config.nx, config.gamma, 1.0 / dx,
                                                               periodic_boundary);
        fp_cuda_add_scaled_1d_kernel<<<grid_size, block_size>>>(d_cells, d_rhs1, d_stage1, config.nx, dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Stage 1 kernels failed in cuda_run_weno2_1d.");
        }

        fp_cuda_weno2_rhs_1d_kernel<<<grid_size, block_size>>>(d_stage1, d_rhs2, config.nx, config.gamma, 1.0 / dx,
                                                               periodic_boundary);
        fp_cuda_rk2_finalize_1d_kernel<<<grid_size, block_size>>>(d_cells, d_stage1, d_rhs2, d_next, config.nx, dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Stage 2 kernels failed in cuda_run_weno2_1d.");
        }

        std::swap(d_cells, d_next);
        time += dt;
        ++step;
    }

    if (cudaMemcpy(result.cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_weno2_1d.");
    }

    cudaFree(d_cells);
    cudaFree(d_rhs1);
    cudaFree(d_stage1);
    cudaFree(d_rhs2);
    cudaFree(d_next);

    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult2D cuda_run_godunov_2d(const RuntimeConfig& config,
                                       const Grid2D& grid,
                                       const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("2D initial state size does not match grid.");
    }

    const std::size_t cell_count = static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    const std::size_t cell_bytes = sizeof(Conservative2D) * cell_count;
    const std::size_t fx_bytes = sizeof(Conservative2D) * static_cast<std::size_t>(grid.nx + 1) *
                                 static_cast<std::size_t>(grid.ny);
    const std::size_t fy_bytes = sizeof(Conservative2D) * static_cast<std::size_t>(grid.nx) *
                                 static_cast<std::size_t>(grid.ny + 1);
    const double hx = dx(grid);
    const double hy = dy(grid);
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative2D* d_cells = nullptr;
    Conservative2D* d_next = nullptr;
    Conservative2D* d_fx = nullptr;
    Conservative2D* d_fy = nullptr;

    if (cudaMalloc(&d_cells, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_fx, fx_bytes) != cudaSuccess ||
        cudaMalloc(&d_fy, fy_bytes) != cudaSuccess) {
        if (d_cells) cudaFree(d_cells);
        if (d_next) cudaFree(d_next);
        if (d_fx) cudaFree(d_fx);
        if (d_fy) cudaFree(d_fy);
        throw std::runtime_error("cudaMalloc failed in cuda_run_godunov_2d.");
    }

    if (cudaMemcpy(d_cells, initial_state.cells.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_cells);
        cudaFree(d_next);
        cudaFree(d_fx);
        cudaFree(d_fy);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_godunov_2d.");
    }

    Field2D host_field = initial_state;
    double time = 0.0;
    int step = 0;

    const dim3 block_flux(16, 16);
    const dim3 grid_flux_x((grid.nx + 1 + block_flux.x - 1) / block_flux.x, (grid.ny + block_flux.y - 1) / block_flux.y);
    const dim3 grid_flux_y((grid.nx + block_flux.x - 1) / block_flux.x, (grid.ny + 1 + block_flux.y - 1) / block_flux.y);
    const dim3 block_cell(16, 16);
    const dim3 grid_cell((grid.nx + block_cell.x - 1) / block_cell.x, (grid.ny + block_cell.y - 1) / block_cell.y);

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemcpy(host_field.cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            throw std::runtime_error("cudaMemcpy D2H state failed in cuda_run_godunov_2d.");
        }

        double dt = fp_host_compute_time_step_2d(host_field, grid, config);
        dt = std::min(dt, config.final_time - time);

        fp_cuda_godunov_flux_x_2d_kernel<<<grid_flux_x, block_flux>>>(d_cells, d_fx, grid.nx, grid.ny, config.gamma,
                                                                      periodic_boundary);
        fp_cuda_godunov_flux_y_2d_kernel<<<grid_flux_y, block_flux>>>(d_cells, d_fy, grid.nx, grid.ny, config.gamma,
                                                                      periodic_boundary);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("2D flux kernels failed in cuda_run_godunov_2d.");
        }

        fp_cuda_godunov_update_2d_kernel<<<grid_cell, block_cell>>>(d_cells, d_fx, d_fy, d_next, grid.nx, grid.ny,
                                                                    dt / hx, dt / hy);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("2D update kernel failed in cuda_run_godunov_2d.");
        }

        std::swap(d_cells, d_next);
        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field.nx = grid.nx;
    result.field.ny = grid.ny;
    result.field.cells.resize(cell_count);
    if (cudaMemcpy(result.field.cells.data(), d_cells, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_godunov_2d.");
    }

    cudaFree(d_cells);
    cudaFree(d_next);
    cudaFree(d_fx);
    cudaFree(d_fy);

    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult2D cuda_run_muscl_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("2D initial state size does not match grid.");
    }

    const std::size_t cell_count = static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    const std::size_t cell_bytes = sizeof(Conservative2D) * cell_count;
    const double inv_dx = 1.0 / dx(grid);
    const double inv_dy = 1.0 / dy(grid);
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative2D* d_field = nullptr;
    Conservative2D* d_rhs1 = nullptr;
    Conservative2D* d_stage1 = nullptr;
    Conservative2D* d_rhs2 = nullptr;
    Conservative2D* d_next = nullptr;

    if (cudaMalloc(&d_field, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_stage1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs2, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess) {
        if (d_field) cudaFree(d_field);
        if (d_rhs1) cudaFree(d_rhs1);
        if (d_stage1) cudaFree(d_stage1);
        if (d_rhs2) cudaFree(d_rhs2);
        if (d_next) cudaFree(d_next);
        throw std::runtime_error("cudaMalloc failed in cuda_run_muscl_2d.");
    }

    if (cudaMemcpy(d_field, initial_state.cells.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_field);
        cudaFree(d_rhs1);
        cudaFree(d_stage1);
        cudaFree(d_rhs2);
        cudaFree(d_next);
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_muscl_2d.");
    }

    Field2D host_field = initial_state;
    double time = 0.0;
    int step = 0;

    const dim3 block2d(16, 16);
    const dim3 grid2d((grid.nx + block2d.x - 1) / block2d.x, (grid.ny + block2d.y - 1) / block2d.y);
    const int block1d = 256;
    const int grid1d = static_cast<int>((cell_count + block1d - 1) / block1d);

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemcpy(host_field.cells.data(), d_field, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
            throw std::runtime_error("cudaMemcpy D2H state failed in cuda_run_muscl_2d.");
        }

        double dt = fp_host_compute_time_step_2d(host_field, grid, config);
        dt = std::min(dt, config.final_time - time);

        fp_cuda_muscl_rhs_2d_kernel<<<grid2d, block2d>>>(d_field, d_rhs1, grid.nx, grid.ny, config.gamma, inv_dx,
                                                         inv_dy, periodic_boundary);
        fp_cuda_add_scaled_2d_kernel<<<grid1d, block1d>>>(d_field, d_rhs1, d_stage1, static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Stage 1 kernels failed in cuda_run_muscl_2d.");
        }

        fp_cuda_muscl_rhs_2d_kernel<<<grid2d, block2d>>>(d_stage1, d_rhs2, grid.nx, grid.ny, config.gamma, inv_dx,
                                                         inv_dy, periodic_boundary);
        fp_cuda_rk2_finalize_2d_kernel<<<grid1d, block1d>>>(d_field, d_stage1, d_rhs2, d_next,
                                                            static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            throw std::runtime_error("Stage 2 kernels failed in cuda_run_muscl_2d.");
        }

        std::swap(d_field, d_next);
        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field.nx = grid.nx;
    result.field.ny = grid.ny;
    result.field.cells.resize(cell_count);
    if (cudaMemcpy(result.field.cells.data(), d_field, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_muscl_2d.");
    }

    cudaFree(d_field);
    cudaFree(d_rhs1);
    cudaFree(d_stage1);
    cudaFree(d_rhs2);
    cudaFree(d_next);

    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult2D cuda_run_weno2_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("2D initial state size does not match grid.");
    }

    const std::size_t cell_count = static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    const std::size_t cell_bytes = sizeof(Conservative2D) * cell_count;
    const double inv_dx = 1.0 / dx(grid);
    const double inv_dy = 1.0 / dy(grid);
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative2D* d_field = nullptr;
    Conservative2D* d_rhs1 = nullptr;
    Conservative2D* d_stage1 = nullptr;
    Conservative2D* d_rhs2 = nullptr;
    Conservative2D* d_next = nullptr;
    double* d_max_x = nullptr;
    double* d_max_y = nullptr;

    auto free_device_memory = [&]() {
        if (d_field) cudaFree(d_field);
        if (d_rhs1) cudaFree(d_rhs1);
        if (d_stage1) cudaFree(d_stage1);
        if (d_rhs2) cudaFree(d_rhs2);
        if (d_next) cudaFree(d_next);
        if (d_max_x) cudaFree(d_max_x);
        if (d_max_y) cudaFree(d_max_y);
    };

    if (cudaMalloc(&d_field, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_stage1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs2, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_max_x, sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_max_y, sizeof(double)) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMalloc failed in cuda_run_weno2_2d.");
    }

    if (cudaMemcpy(d_field, initial_state.cells.data(), cell_bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_weno2_2d.");
    }

    double time = 0.0;
    int step = 0;

    const dim3 block2d(16, 16);
    const dim3 grid2d((grid.nx + block2d.x - 1) / block2d.x, (grid.ny + block2d.y - 1) / block2d.y);
    const int block1d = 256;
    const int grid1d = static_cast<int>((cell_count + block1d - 1) / block1d);

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemset(d_max_x, 0, sizeof(double)) != cudaSuccess ||
            cudaMemset(d_max_y, 0, sizeof(double)) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("CFL reduction reset failed in cuda_run_weno2_2d.");
        }

        fp_cuda_max_signal_speed_2d_kernel<<<grid1d, block1d>>>(
            d_field, static_cast<int>(cell_count), config.gamma, d_max_x, d_max_y);

        double max_x = 0.0;
        double max_y = 0.0;
        if (cudaMemcpy(&max_x, d_max_x, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess ||
            cudaMemcpy(&max_y, d_max_y, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("CFL reduction copy failed in cuda_run_weno2_2d.");
        }

        const double max_speed =
            std::max(std::max(max_x, max_y), 1.0e-12);
        double dt =
            config.cfl * std::min(dx(grid), dy(grid)) / max_speed;
        dt = std::min(dt, config.final_time - time);

        fp_cuda_weno2_rhs_2d_kernel<<<grid2d, block2d>>>(d_field, d_rhs1, grid.nx, grid.ny, config.gamma, inv_dx,
                                                         inv_dy, periodic_boundary);
        fp_cuda_add_scaled_2d_kernel<<<grid1d, block1d>>>(d_field, d_rhs1, d_stage1, static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("Stage 1 kernels failed in cuda_run_weno2_2d.");
        }

        fp_cuda_weno2_rhs_2d_kernel<<<grid2d, block2d>>>(d_stage1, d_rhs2, grid.nx, grid.ny, config.gamma, inv_dx,
                                                         inv_dy, periodic_boundary);
        fp_cuda_rk2_finalize_2d_kernel<<<grid1d, block1d>>>(d_field, d_stage1, d_rhs2, d_next,
                                                            static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("Stage 2 kernels failed in cuda_run_weno2_2d.");
        }

        std::swap(d_field, d_next);
        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field.nx = grid.nx;
    result.field.ny = grid.ny;
    result.field.cells.resize(cell_count);
    if (cudaMemcpy(result.field.cells.data(), d_field, cell_bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_weno2_2d.");
    }

    free_device_memory();

    result.final_time = time;
    result.steps = step;
    return result;
}

SimulationResult2D cuda_run_weno3_2d(const RuntimeConfig& config,
                                     const Grid2D& grid,
                                     const Field2D& initial_state) {
    if (initial_state.nx != grid.nx || initial_state.ny != grid.ny) {
        throw std::runtime_error("2D initial state size does not match grid.");
    }

    const std::size_t cell_count =
        static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    const std::size_t cell_bytes = sizeof(Conservative2D) * cell_count;
    const double inv_dx = 1.0 / dx(grid);
    const double inv_dy = 1.0 / dy(grid);
    const int periodic_boundary = (config.boundary == "periodic") ? 1 : 0;

    Conservative2D* d_field = nullptr;
    Conservative2D* d_rhs1 = nullptr;
    Conservative2D* d_stage1 = nullptr;
    Conservative2D* d_rhs2 = nullptr;
    Conservative2D* d_stage2 = nullptr;
    Conservative2D* d_rhs3 = nullptr;
    Conservative2D* d_next = nullptr;
    double* d_max_x = nullptr;
    double* d_max_y = nullptr;

    auto free_device_memory = [&]() {
        if (d_field) cudaFree(d_field);
        if (d_rhs1) cudaFree(d_rhs1);
        if (d_stage1) cudaFree(d_stage1);
        if (d_rhs2) cudaFree(d_rhs2);
        if (d_stage2) cudaFree(d_stage2);
        if (d_rhs3) cudaFree(d_rhs3);
        if (d_next) cudaFree(d_next);
        if (d_max_x) cudaFree(d_max_x);
        if (d_max_y) cudaFree(d_max_y);
    };

    if (cudaMalloc(&d_field, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_stage1, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs2, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_stage2, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_rhs3, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_next, cell_bytes) != cudaSuccess ||
        cudaMalloc(&d_max_x, sizeof(double)) != cudaSuccess ||
        cudaMalloc(&d_max_y, sizeof(double)) != cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMalloc failed in cuda_run_weno3_2d.");
    }

    if (cudaMemcpy(d_field, initial_state.cells.data(), cell_bytes, cudaMemcpyHostToDevice) !=
        cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("cudaMemcpy H2D failed in cuda_run_weno3_2d.");
    }

    double time = 0.0;
    int step = 0;

    const dim3 block2d(16, 16);
    const dim3 grid2d(
        (grid.nx + block2d.x - 1) / block2d.x,
        (grid.ny + block2d.y - 1) / block2d.y);
    const int block1d = 256;
    const int grid1d = static_cast<int>((cell_count + block1d - 1) / block1d);

    while (time < config.final_time && step < config.max_steps) {
        if (cudaMemset(d_max_x, 0, sizeof(double)) != cudaSuccess ||
            cudaMemset(d_max_y, 0, sizeof(double)) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("CFL reduction reset failed in cuda_run_weno3_2d.");
        }

        fp_cuda_max_signal_speed_2d_kernel<<<grid1d, block1d>>>(
            d_field, static_cast<int>(cell_count), config.gamma, d_max_x, d_max_y);

        double max_x = 0.0;
        double max_y = 0.0;
        if (cudaMemcpy(&max_x, d_max_x, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess ||
            cudaMemcpy(&max_y, d_max_y, sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("CFL reduction copy failed in cuda_run_weno3_2d.");
        }

        const double spectral_sum =
            std::max(max_x, 1.0e-12) * inv_dx + std::max(max_y, 1.0e-12) * inv_dy;
        double dt = config.cfl / spectral_sum;
        dt = std::min(dt, config.final_time - time);

        fp_cuda_weno3_rhs_2d_kernel<<<grid2d, block2d>>>(
            d_field, d_rhs1, grid.nx, grid.ny, config.gamma, inv_dx, inv_dy,
            periodic_boundary);
        fp_cuda_add_scaled_2d_kernel<<<grid1d, block1d>>>(
            d_field, d_rhs1, d_stage1, static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("Stage 1 kernels failed in cuda_run_weno3_2d.");
        }

        fp_cuda_weno3_rhs_2d_kernel<<<grid2d, block2d>>>(
            d_stage1, d_rhs2, grid.nx, grid.ny, config.gamma, inv_dx, inv_dy,
            periodic_boundary);
        fp_cuda_rk3_stage2_2d_kernel<<<grid1d, block1d>>>(
            d_field, d_stage1, d_rhs2, d_stage2, static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("Stage 2 kernels failed in cuda_run_weno3_2d.");
        }

        fp_cuda_weno3_rhs_2d_kernel<<<grid2d, block2d>>>(
            d_stage2, d_rhs3, grid.nx, grid.ny, config.gamma, inv_dx, inv_dy,
            periodic_boundary);
        fp_cuda_rk3_finalize_2d_kernel<<<grid1d, block1d>>>(
            d_field, d_stage2, d_rhs3, d_next, static_cast<int>(cell_count), dt);
        if (cudaDeviceSynchronize() != cudaSuccess) {
            free_device_memory();
            throw std::runtime_error("Stage 3 kernels failed in cuda_run_weno3_2d.");
        }

        std::swap(d_field, d_next);
        time += dt;
        ++step;
    }

    SimulationResult2D result;
    result.grid = grid;
    result.field.nx = grid.nx;
    result.field.ny = grid.ny;
    result.field.cells.resize(cell_count);
    if (cudaMemcpy(result.field.cells.data(), d_field, cell_bytes, cudaMemcpyDeviceToHost) !=
        cudaSuccess) {
        free_device_memory();
        throw std::runtime_error("Final cudaMemcpy D2H failed in cuda_run_weno3_2d.");
    }

    free_device_memory();
    result.final_time = time;
    result.steps = step;
    return result;
}

}  // namespace fp
