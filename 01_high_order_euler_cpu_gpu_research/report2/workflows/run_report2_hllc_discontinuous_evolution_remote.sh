#!/usr/bin/env bash
set -euo pipefail

package_name="report2_hllc_discontinuous_evolution_source.tar.gz"
work_dir="${HOME}/report2_hllc_discontinuous_evolution_20260725"
result_archive="${HOME}/report2_hllc_discontinuous_evolution_results.tar.gz"
resolution="${REPORT2_HLLC_EVOLUTION_RESOLUTION:-600}"
result_root="results/report2_hllc/discontinuous_evolution_2d"
output_dir="${result_root}/fields"

rm -rf "${work_dir}"
rm -f "${result_archive}"
mkdir -p "${work_dir}"
tar -xzf "${HOME}/${package_name}" -C "${work_dir}"
cd "${work_dir}"

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export LD_LIBRARY_PATH="/lsc/opt/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export TMPDIR="${work_dir}/.compiler_tmp"
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
mkdir -p "${TMPDIR}" "${output_dir}"

echo "[1/5] Configuring the HLLC CUDA build"
cmake -S . -B build-report2-hllc-evolution \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release

echo "[2/5] Building the solver"
cmake --build build-report2-hllc-evolution -j 8

binary="./build-report2-hllc-evolution/final_project"

echo "[3/5] Mandatory CUDA smoke test"
"${binary}" \
    --backend=cuda \
    --problem=quadrant_riemann_2d \
    --method=weno3 \
    --nx=25 \
    --ny=25 \
    --boundary=outflow \
    --final_time=0.002 \
    --cfl=0.2 \
    --output=none

run_snapshot() {
    local problem="$1"
    local final_time="$2"
    local time_tag="$3"
    local output="${output_dir}/${problem}_weno3_cuda_n${resolution}_${time_tag}.csv"

    echo "  problem=${problem} method=weno3 n=${resolution} t=${final_time}"
    "${binary}" \
        --backend=cuda \
        --problem="${problem}" \
        --method=weno3 \
        --nx="${resolution}" \
        --ny="${resolution}" \
        --boundary=outflow \
        --final_time="${final_time}" \
        --cfl=0.2 \
        --output="${output}"
}

echo "[4/5] Exporting initial and intermediate discontinuous fields"
for problem in quadrant_riemann_2d shock_bubble_2d; do
    run_snapshot "${problem}" 0.0 t000
    run_snapshot "${problem}" 0.1 t010
done

echo "[5/5] Verifying and packaging snapshots"
REPORT2_HLLC_EVOLUTION_RESOLUTION="${resolution}" python3 - <<'PY'
import csv
import math
import os
from pathlib import Path

root = Path("results/report2_hllc/discontinuous_evolution_2d/fields")
paths = sorted(root.glob("*.csv"))
resolution = int(os.environ["REPORT2_HLLC_EVOLUTION_RESOLUTION"])
if len(paths) != 4:
    raise SystemExit(f"expected 4 snapshot CSV files, found {len(paths)}")

for path in paths:
    rows = 0
    min_rho = math.inf
    min_p = math.inf
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            values = [float(row[name]) for name in ("rho", "u", "v", "p", "energy")]
            if not all(math.isfinite(value) for value in values):
                raise SystemExit(f"non-finite state in {path}")
            min_rho = min(min_rho, float(row["rho"]))
            min_p = min(min_p, float(row["p"]))
            rows += 1
    if rows != resolution * resolution:
        raise SystemExit(f"unexpected row count in {path}: {rows}")
    if min_rho <= 0.0 or min_p <= 0.0:
        raise SystemExit(
            f"non-positive state in {path}: min_rho={min_rho}, min_p={min_p}"
        )
    print(
        f"  PASS {path.name}: rows={rows}, "
        f"min_rho={min_rho:.6e}, min_p={min_p:.6e}"
    )
PY

tar -czf "${result_archive}" -C "${work_dir}" "${result_root}"
echo "result archive: ${result_archive}"
