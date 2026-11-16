#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_weno2_openmp_refresh_source.tar.gz"
work_dir="${HOME}/report2_hllc_weno2_openmp_refresh_20260717"
result_archive="${HOME}/report2_hllc_weno2_openmp_refresh_results.tar.gz"
result_dir="results/report2_hllc/performance/weno2_openmp_refresh"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

mkdir -p "${work_dir}/.compiler_tmp"
export TMPDIR="${work_dir}/.compiler_tmp"
export OMP_NUM_THREADS=16
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${result_dir}"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "OMP_DYNAMIC=${OMP_DYNAMIC}"
    echo "OMP_PROC_BIND=${OMP_PROC_BIND}"
    echo "OMP_PLACES=${OMP_PLACES}"
    lscpu | sed -n '1,30p'
} > "${result_dir}/environment.txt" 2>&1

echo "[1/7] Configuring Release OpenMP build"
cmake -S . -B build-report2-hllc-openmp \
    -DFP_ENABLE_CUDA=OFF \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release

echo "[2/7] Building solver"
cmake --build build-report2-hllc-openmp -j 8

echo "[3/7] Refreshing smooth WENO-2 OpenMP scaling"
python3 scripts/run_report2_runtime_scaling_2d.py \
    --build-dir build-report2-hllc-openmp \
    --methods weno2 \
    --backends openmp \
    --resolutions 25 50 100 200 600 \
    --warmup 1 \
    --repeats 3 \
    --openmp-threads 16 \
    --output "${result_dir}/smooth_weno2_openmp_refresh.csv"

echo "[4/7] Refreshing discontinuous WENO-2 OpenMP n=50,100,200"
python3 scripts/run_report2_discontinuous_runtime_2d.py \
    --build-dir build-report2-hllc-openmp \
    --methods weno2 \
    --backends openmp \
    --resolutions 50 100 200 \
    --warmup 1 \
    --repeats 3 \
    --openmp-threads 16 \
    --output "${result_dir}/discontinuous_weno2_openmp_refresh.csv"

echo "[5/7] Refreshing discontinuous WENO-2 OpenMP n=600"
python3 scripts/run_report2_discontinuous_runtime_2d.py \
    --build-dir build-report2-hllc-openmp \
    --methods weno2 \
    --backends openmp \
    --resolutions 600 \
    --warmup 0 \
    --repeats 1 \
    --openmp-threads 16 \
    --output "${result_dir}/discontinuous_weno2_openmp_refresh.csv" \
    --resume

echo "[6/7] Refreshing WENO-2 OpenMP thread scaling and phase profiles"
python3 scripts/run_report2_openmp_thread_scaling_2d.py \
    --build-dir build-report2-hllc-openmp \
    --methods weno2 \
    --threads 1 2 4 8 16 \
    --nx 100 \
    --ny 100 \
    --warmup 1 \
    --repeats 3 \
    --output "${result_dir}/openmp_thread_scaling_weno2_refresh.csv"

python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-openmp \
    --resolutions 100 \
    --methods weno2 \
    --backends openmp \
    --repeats 3 \
    --openmp-threads 16 \
    --output "${result_dir}/phase_profile_weno2_openmp_refresh.csv"

python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-openmp \
    --resolutions 600 \
    --methods weno2 \
    --backends openmp \
    --repeats 1 \
    --openmp-threads 16 \
    --output "${result_dir}/phase_profile_weno2_openmp_refresh.csv" \
    --resume

echo "[7/7] Packaging refresh candidates"
tar -czf "${result_archive}" -C "${work_dir}" "${result_dir}"
echo "result archive: ${result_archive}"
