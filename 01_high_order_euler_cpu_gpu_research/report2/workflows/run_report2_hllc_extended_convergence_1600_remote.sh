#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_extended_convergence_1600_source.tar.gz"
work_dir="${HOME}/report2_hllc_extended_convergence_1600_20260721"
result_archive="${HOME}/report2_hllc_extended_convergence_1600_results.tar.gz"
result_root="results/report2_hllc/extended_convergence_1600"
build_dir="build-report2-hllc-extended"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export TMPDIR="${work_dir}/.compiler_tmp"
mkdir -p "${TMPDIR}" "${result_root}/validation_1d_openmp" \
    "${result_root}/validation_2d_cuda"

echo "[1/6] Recording environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    /lsc/opt/cuda-12.9/bin/nvcc --version
} > "${result_root}/environment.txt" 2>&1

echo "[2/6] Configuring and building"
cmake -S . -B "${build_dir}" \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "${build_dir}" -j 8

echo "[3/6] Mandatory CUDA smoke test"
"${build_dir}/final_project" \
    --backend=cuda --problem=gaussian_advection_2d --method=weno2 \
    --nx=25 --ny=25 --boundary=periodic --final_time=0.002 \
    --cfl=0.4 --output=none

echo "[4/6] Extending 1D OpenMP convergence to n=1600"
python3 scripts/run_report2_validation_1d.py \
    --build-dir "${build_dir}" \
    --output-dir "${result_root}/validation_1d_openmp" \
    --backend openmp --resolutions 800 1600 --final-time 0.3

echo "[5/6] Extending 2D CUDA convergence to n=1600"
python3 scripts/run_report2_extended_convergence_2d.py \
    --build-dir "${build_dir}" \
    --output-dir "${result_root}/validation_2d_cuda" \
    --resolutions 400 800 1600 --final-time 0.1

echo "[6/6] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" "${result_root}"
echo "result archive: ${result_archive}"
