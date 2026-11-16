#!/usr/bin/env bash
set -euo pipefail

work_dir="${HOME}/report2_hllc_phase_profile_20260714"
result_archive="${HOME}/report2_hllc_phase_profile_results.tar.gz"
result_dir="results/report2_hllc/performance/phase_profile"
result_csv="${result_dir}/phase_profile_2d_hllc_lovelace_a30.csv"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export LD_LIBRARY_PATH="/lsc/opt/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export TMPDIR="${work_dir}/.compiler_tmp"
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OMP_PLACES=cores

if [[ ! -f "${work_dir}/${result_csv}" ]]; then
    echo "Existing partial CSV was not found: ${work_dir}/${result_csv}" >&2
    exit 2
fi

cd "${work_dir}"
mkdir -p "${TMPDIR}"
rm -f "${result_archive}"

echo "[resume 1/5] Reconfiguring with Lovelace CUDA 12.9"
echo "nvcc: $(command -v nvcc)"
nvcc --version
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release

echo "[resume 2/5] Rebuilding solver"
cmake --build build-report2-hllc-cuda -j 8

echo "[resume 3/5] Completing 100 x 100 profiles"
python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolutions 100 \
    --methods weno2 weno3 ader3 \
    --backends serial openmp cuda \
    --repeats 3 \
    --openmp-threads 16 \
    --output "${result_csv}" \
    --resume

echo "[resume 4/5] Completing 600 x 600 profiles"
python3 scripts/run_report2_phase_profile_2d.py \
    --build-dir build-report2-hllc-cuda \
    --resolutions 600 \
    --methods weno2 weno3 ader3 \
    --backends serial openmp cuda \
    --repeats 1 \
    --openmp-threads 16 \
    --output "${result_csv}" \
    --resume

echo "[resume 5/5] Plotting, summarizing, and packaging"
python3 scripts/plot_report2_phase_profile_2d.py "${result_csv}"
python3 scripts/summarize_report2_phase_profile_2d.py "${result_csv}"
tar -czf "${result_archive}" -C "${work_dir}" "${result_dir}"
echo "result archive: ${result_archive}"
