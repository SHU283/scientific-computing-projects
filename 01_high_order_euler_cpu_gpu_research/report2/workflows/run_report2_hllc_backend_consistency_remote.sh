#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_backend_consistency_source.tar.gz"
work_dir="${HOME}/report2_hllc_backend_consistency_20260625"
result_archive="${HOME}/report2_hllc_backend_consistency_results.tar.gz"
result_dir="results/report2_hllc/backend_consistency_2d_lovelace"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${result_dir}"

echo "[1/6] Recording Lovelace environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "nvcc: $(command -v nvcc)"
    nvcc --version
    nvidia-smi || echo "warning: nvidia-smi returned a non-zero status"
} > "${result_dir}/environment.txt" 2>&1

echo "[2/6] Configuring CUDA/OpenMP build"
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES=80

echo "[3/6] Building solver"
cmake --build build-report2-hllc-cuda -j 8

echo "[4/6] Running representative HLLC backend-consistency checks"
python3 scripts/run_report2_backend_consistency_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolution 100 \
    --openmp-threads "${OMP_NUM_THREADS}" \
    --output-dir "${result_dir}"

echo "[5/6] Plotting backend-consistency summary"
python3 scripts/plot_report2_backend_consistency_2d.py \
    "${result_dir}/backend_consistency_2d_report2_lovelace_a30.csv"

echo "[6/6] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" "${result_dir}"
echo "result archive: ${result_archive}"
