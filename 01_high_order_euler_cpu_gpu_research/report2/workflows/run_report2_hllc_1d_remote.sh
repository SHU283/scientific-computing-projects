#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_1d_source.tar.gz"
work_dir="${HOME}/report2_hllc_1d_20260702"
result_archive="${HOME}/report2_hllc_1d_results.tar.gz"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

echo "[1/5] Configuring the Report 2 HLLC CPU/OpenMP build"
cmake -S . -B build-report2-hllc \
    -DFP_ENABLE_CUDA=OFF \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release

echo "[2/5] Building the solver"
cmake --build build-report2-hllc -j 8

echo "[3/5] Running 1D smooth convergence validation"
python3 scripts/run_report2_validation_1d.py \
    --build-dir build-report2-hllc \
    --output-dir results/report2_hllc/validation_1d \
    --resolutions 50 100 200 400 800 \
    --final-time 0.3 \
    --backend openmp

echo "[4/5] Running highest-resolution 1D field exports"
python3 scripts/run_report2_high_resolution_fields_1d.py \
    --build-dir build-report2-hllc \
    --output-dir results/report2_hllc/high_resolution_fields/1d \
    --resolution 800

echo "[5/5] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" \
    results/report2_hllc/validation_1d \
    results/report2_hllc/high_resolution_fields/1d

echo "result archive: ${result_archive}"
