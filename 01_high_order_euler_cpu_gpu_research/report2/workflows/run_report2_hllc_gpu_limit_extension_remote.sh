#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_gpu_limit_extension_source.tar.gz"
work_dir="${HOME}/report2_hllc_gpu_limit_extension_20260702"
result_archive="${HOME}/report2_hllc_gpu_limit_extension_results.tar.gz"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_DYNAMIC=FALSE

output_dir="results/report2_hllc/performance/gpu_limit"
runtime_csv="${output_dir}/runtime_accuracy_2d_report2_hllc_cuda_gpu_limit_extension_lovelace_a30.csv"
mkdir -p "${output_dir}"

echo "[1/4] Recording environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}"
    echo "nvcc: $(command -v nvcc)"
    nvcc --version
    nvidia-smi || echo "warning: nvidia-smi returned a non-zero status"
} > "${output_dir}/environment_gpu_limit_extension.txt" 2>&1

echo "[2/4] Configuring and building CUDA solver"
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES=80
cmake --build build-report2-hllc-cuda -j 8

echo "[3/4] Running HLLC CUDA-only high-resolution extension sweep"
python3 scripts/run_report2_gpu_limit_extension_2d.py \
    --build-dir build-report2-hllc-cuda \
    --methods weno2 weno3 ader3 \
    --resolutions 2400 3000 3600 4200 4600 5000 \
    --warmup 0 \
    --repeats 1 \
    --output "${runtime_csv}" \
    --failure-output "${output_dir}/runtime_accuracy_2d_report2_hllc_cuda_gpu_limit_extension_failures_lovelace_a30.csv" \
    --keep-going \
    --resume

echo "[4/4] Packaging incremental results"
tar -czf "${result_archive}" -C "${work_dir}" "${output_dir}"
echo "result archive: ${result_archive}"
