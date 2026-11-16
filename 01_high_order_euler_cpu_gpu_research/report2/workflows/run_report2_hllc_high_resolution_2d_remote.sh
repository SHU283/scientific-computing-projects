#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_high_resolution_2d_source.tar.gz"
work_dir="${HOME}/report2_hllc_high_resolution_2d_20260628"
result_archive="${HOME}/report2_hllc_high_resolution_2d_results.tar.gz"
resolution="${REPORT2_HLLC_FIELD_RESOLUTION:-600}"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"

output_dir="results/report2_hllc/high_resolution_fields/2d"
mkdir -p "${output_dir}"

echo "[1/5] Configuring the Report 2 HLLC CUDA build"
cmake -S . -B build-report2-hllc-cuda \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CUDA_ARCHITECTURES=80

echo "[2/5] Building the solver"
cmake --build build-report2-hllc-cuda -j 8

run_field() {
    local problem="$1"
    local method="$2"
    local boundary="$3"
    local final_time="$4"
    local cfl="$5"
    local output="${output_dir}/${problem}_${method}_n${resolution}.csv"

    echo "  problem=${problem} method=${method} n=${resolution} t=${final_time} cfl=${cfl}"
    ./build-report2-hllc-cuda/final_project \
        --backend=cuda \
        --problem="${problem}" \
        --method="${method}" \
        --nx="${resolution}" \
        --ny="${resolution}" \
        --boundary="${boundary}" \
        --final_time="${final_time}" \
        --cfl="${cfl}" \
        --output="${output}"
}

echo "[3/5] Running the smooth 2D Gaussian fields"
run_field gaussian_advection_2d weno2 periodic 0.1 0.4
run_field gaussian_advection_2d weno3 periodic 0.1 0.4
run_field gaussian_advection_2d ader3 periodic 0.1 0.25

echo "[4/5] Running the discontinuous 2D fields"
for method in weno2 weno3; do
    run_field quadrant_riemann_2d "${method}" outflow 0.2 0.2
    run_field shock_bubble_2d "${method}" outflow 0.2 0.2
done
run_field quadrant_riemann_2d ader3 outflow 0.2 0.14
run_field shock_bubble_2d ader3 outflow 0.2 0.14

echo "[5/5] Checking exported fields and packaging results"
python3 - <<'PY'
import csv
import math
from pathlib import Path

root = Path("results/report2_hllc/high_resolution_fields/2d")
paths = sorted(root.glob("*.csv"))
if len(paths) != 9:
    raise SystemExit(f"expected 9 field CSV files, found {len(paths)}")

for path in paths:
    row_count = 0
    min_rho = math.inf
    min_p = math.inf
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            values = [float(row[name]) for name in ("rho", "u", "v", "p", "energy")]
            if not all(math.isfinite(value) for value in values):
                raise SystemExit(f"non-finite value in {path}")
            min_rho = min(min_rho, float(row["rho"]))
            min_p = min(min_p, float(row["p"]))
            row_count += 1
    if min_rho <= 0.0 or min_p <= 0.0:
        raise SystemExit(
            f"non-positive state in {path}: min_rho={min_rho}, min_p={min_p}"
        )
    print(
        f"  PASS {path.name}: rows={row_count}, "
        f"min_rho={min_rho:.6e}, min_p={min_p:.6e}"
    )
PY

tar -czf "${result_archive}" -C "${work_dir}" "${output_dir}"
echo "result archive: ${result_archive}"
