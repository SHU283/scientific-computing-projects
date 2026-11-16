#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_convergence_2d_source.tar.gz"
work_dir="${HOME}/report2_hllc_convergence_2d_20260625"
result_archive="${HOME}/report2_hllc_convergence_2d_results.tar.gz"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

output_dir="results/report2_hllc/validation_2d_lovelace"
mkdir -p "${output_dir}"

echo "[1/5] Recording Lovelace environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    lscpu | sed -n '1,20p'
} > "${output_dir}/environment.txt" 2>&1

echo "[2/5] Configuring HLLC CPU/OpenMP build"
cmake -S . -B build-hllc-report2 \
    -DFP_ENABLE_CUDA=OFF \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release

echo "[3/5] Building solver"
cmake --build build-hllc-report2 -j 8

echo "[4/5] Running HLLC 2D Gaussian convergence"
python3 scripts/run_report2_validation_2d.py \
    --build-dir build-hllc-report2 \
    --output-dir "${output_dir}" \
    --methods weno2 weno3 ader3 \
    --backend cuda \
    --resolutions 25 50 100 200 400 \
    --final-time 0.1

python3 scripts/plot_report2_convergence_2d.py \
    "${output_dir}/report2_high_order_convergence_2d.csv" \
    --metric l1_rho \
    --output "${output_dir}/report2_hllc_high_order_convergence_2d.l1_rho.png"

echo "[5/5] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" "${output_dir}"
echo "result archive: ${result_archive}"
