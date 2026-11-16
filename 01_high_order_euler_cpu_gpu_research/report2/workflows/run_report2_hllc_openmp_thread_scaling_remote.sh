#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_openmp_thread_scaling_source.tar.gz"
work_dir="${HOME}/report2_hllc_openmp_thread_scaling_20260702"
result_archive="${HOME}/report2_hllc_openmp_thread_scaling_results.tar.gz"
result_dir="results/report2_hllc/performance/openmp_thread_scaling"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${result_dir}"

echo "[1/5] Recording Lovelace environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "OMP_DYNAMIC=${OMP_DYNAMIC}"
    echo "OMP_PROC_BIND=${OMP_PROC_BIND}"
    echo "OMP_PLACES=${OMP_PLACES}"
    lscpu | sed -n '1,25p'
} > "${result_dir}/environment_openmp_thread_scaling.txt" 2>&1

echo "[2/5] Configuring HLLC OpenMP build"
cmake -S . -B build-report2-hllc-openmp \
    -DFP_ENABLE_CUDA=OFF \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release

echo "[3/5] Building solver"
cmake --build build-report2-hllc-openmp -j 8

echo "[4/5] Running HLLC OpenMP thread scaling"
python3 scripts/run_report2_openmp_thread_scaling_2d.py \
    --build-dir build-report2-hllc-openmp \
    --nx 100 \
    --ny 100 \
    --warmup 1 \
    --repeats 3 \
    --threads 1 2 4 8 16 \
    --methods weno2 weno3 ader3 \
    --output "${result_dir}/openmp_thread_scaling_report2_hllc_lovelace_a30.csv" \
    --resume

python3 scripts/plot_openmp_thread_scaling.py \
    "${result_dir}/openmp_thread_scaling_report2_hllc_lovelace_a30.csv" \
    --metric best_seconds \
    --figsize 10.8 4.2

echo "[5/5] Packaging results"
tar -czf "${result_archive}" -C "${work_dir}" "${result_dir}"
echo "result archive: ${result_archive}"
