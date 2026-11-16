#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_performance_source.tar.gz"
work_dir="${HOME}/report2_hllc_performance_20260628"
result_archive="${HOME}/report2_hllc_performance_results.tar.gz"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

output_dir="results/report2_hllc/performance"
runtime_csv="${output_dir}/runtime_accuracy_2d_report2_hllc_lovelace_a30.csv"
mkdir -p "${output_dir}"

echo "[1/5] Recording environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
    echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}"
    echo "nvcc: $(command -v nvcc)"
    nvcc --version
    nvidia-smi || echo "warning: nvidia-smi returned a non-zero status"
} > "${output_dir}/environment.txt" 2>&1

echo "[2/5] Configuring the CUDA/OpenMP build"
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES=80

echo "[3/5] Building the solver"
cmake --build build-report2-hllc-cuda -j 8

echo "[4/5] Running HLLC smooth 2D no-I/O performance campaign"
python3 scripts/run_report2_runtime_scaling_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolutions 25 50 100 200 600 \
    --warmup 1 \
    --repeats 3 \
    --openmp-threads "${OMP_NUM_THREADS}" \
    --output "${runtime_csv}"

echo "[5/5] Plotting and packaging"
python3 scripts/plot_report2_runtime_scaling_2d.py "${runtime_csv}"
python3 scripts/plot_report2_accuracy_runtime_2d.py "${runtime_csv}"

tar -czf "${result_archive}" -C "${work_dir}" "${output_dir}"
echo "result archive: ${result_archive}"
