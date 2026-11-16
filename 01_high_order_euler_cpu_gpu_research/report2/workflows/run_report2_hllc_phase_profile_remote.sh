#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_phase_profile_source.tar.gz"
work_dir="${HOME}/report2_hllc_phase_profile_20260714"
result_archive="${HOME}/report2_hllc_phase_profile_results.tar.gz"
result_dir="results/report2_hllc/performance/phase_profile"
result_csv="${result_dir}/phase_profile_2d_hllc_lovelace_a30.csv"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export LD_LIBRARY_PATH="/lsc/opt/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

mkdir -p "${result_dir}"

echo "[1/6] Recording Lovelace environment"
{
    echo "hostname: $(hostname)"
    echo "date: $(date --iso-8601=seconds)"
    echo "OMP_DYNAMIC=${OMP_DYNAMIC}"
    echo "OMP_PROC_BIND=${OMP_PROC_BIND}"
    echo "OMP_PLACES=${OMP_PLACES}"
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}"
    echo "nvcc: $(command -v nvcc)"
    nvcc --version
    lscpu | sed -n '1,30p'
    nvidia-smi || true
} > "${result_dir}/environment_phase_profile.txt" 2>&1

echo "[2/6] Configuring HLLC CUDA/OpenMP build"
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release

echo "[3/6] Building solver"
cmake --build build-report2-hllc-cuda -j 8

echo "[4/6] Profiling the low-resolution 100 x 100 case (three repeats)"
python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolutions 100 \
    --methods weno2 weno3 ader3 \
    --backends serial openmp cuda \
    --repeats 3 \
    --openmp-threads 16 \
    --output "${result_csv}"

echo "[5/6] Profiling the high-resolution 600 x 600 case (one repeat)"
python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolutions 600 \
    --methods weno2 weno3 ader3 \
    --backends serial openmp cuda \
    --repeats 1 \
    --openmp-threads 16 \
    --output "${result_csv}" \
    --resume

echo "[6/6] Plotting and packaging phase profiles"
python3 scripts/plot_report2_phase_profile_2d.py "${result_csv}"
python3 scripts/summarize_report2_phase_profile_2d.py "${result_csv}"
tar -czf "${result_archive}" -C "${work_dir}" "${result_dir}"
echo "result archive: ${result_archive}"
