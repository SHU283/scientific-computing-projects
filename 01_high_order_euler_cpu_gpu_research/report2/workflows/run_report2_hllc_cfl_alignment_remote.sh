#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_cfl_alignment_source.tar.gz"
work_dir="${HOME}/report2_hllc_cfl_alignment_20260720"
result_archive="${HOME}/report2_hllc_cfl_alignment_results.tar.gz"
result_root="results/report2_hllc/cfl_aligned_validation"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export TMPDIR="${work_dir}/.compiler_tmp"
mkdir -p "${TMPDIR}" "${result_root}/validation_2d" \
    "${result_root}/validation_1d" "${result_root}/high_resolution_fields_1d" \
    "${result_root}/backend_consistency_2d_lovelace"

echo "[1/8] Recording environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    /lsc/opt/cuda-12.9/bin/nvcc --version
} > "${result_root}/environment.txt" 2>&1

echo "[2/8] Configuring and building HLLC solver"
cmake -S . -B build-report2-hllc-cfl-audit \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build-report2-hllc-cfl-audit -j 8

echo "[3/8] Mandatory CUDA smoke test"
build-report2-hllc-cfl-audit/final_project \
    --backend=cuda --problem=gaussian_advection_2d --method=weno2 \
    --nx=25 --ny=25 --boundary=periodic --final_time=0.002 \
    --cfl=0.4 --output=none

echo "[4/8] Official 2D smooth convergence"
python3 scripts/run_report2_validation_2d.py \
    --build-dir build-report2-hllc-cfl-audit \
    --output-dir "${result_root}/validation_2d" \
    --backend cuda --methods weno2 weno3 ader3 \
    --resolutions 25 50 100 200 400 --final-time 0.1
python3 scripts/plot_report2_convergence_2d.py \
    "${result_root}/validation_2d/report2_high_order_convergence_2d.csv" \
    --metric l1_rho \
    --output "${result_root}/validation_2d/report2_high_order_convergence_2d.l1_rho.png"

echo "[5/8] Official 1D smooth convergence"
python3 scripts/run_report2_validation_1d.py \
    --build-dir build-report2-hllc-cfl-audit \
    --output-dir "${result_root}/validation_1d" \
    --backend openmp --resolutions 50 100 200 400 800 --final-time 0.3
python3 scripts/plot_report2_convergence_1d.py \
    "${result_root}/validation_1d/report2_high_order_convergence_1d.csv" \
    --output "${result_root}/validation_1d/report2_high_order_convergence_1d.png"

echo "[6/8] CFL-aligned 1D field exports"
python3 scripts/run_report2_high_resolution_fields_1d.py \
    --build-dir build-report2-hllc-cfl-audit \
    --output-dir "${result_root}/high_resolution_fields_1d" \
    --resolution 800

echo "[7/8] CFL-aligned backend consistency"
python3 scripts/run_report2_backend_consistency_2d.py \
    --build-dir build-report2-hllc-cfl-audit \
    --resolution 100 --openmp-threads "${OMP_NUM_THREADS}" \
    --output-dir "${result_root}/backend_consistency_2d_lovelace"
python3 scripts/plot_report2_backend_consistency_2d.py \
    "${result_root}/backend_consistency_2d_lovelace/backend_consistency_2d_report2_lovelace_a30.csv"

echo "[8/8] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" "${result_root}"
echo "result archive: ${result_archive}"
