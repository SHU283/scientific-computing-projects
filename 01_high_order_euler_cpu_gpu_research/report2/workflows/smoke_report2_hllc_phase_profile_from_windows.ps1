param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$Root = "D:\code\Final_projects"
$ProbeSource = Join-Path $Root "report2-hllc\src\backends\cuda_backend.cu"
$RemoteProbeSource = "~/report2_hllc_phase_profile_20260714/src/backends/cuda_backend.cu"

Write-Host "=== Uploading diagnostic CUDA probe source ==="
scp $ProbeSource "${User}@${HostName}:$RemoteProbeSource"
if ($LASTEXITCODE -ne 0) {
    throw "CUDA probe source upload failed with exit code $LASTEXITCODE"
}

$remote = @'
set -euo pipefail

export PATH="/lsc/opt/cuda-12.9/bin:${PATH}"
export LD_LIBRARY_PATH="/lsc/opt/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES=0
export FP_PHASE_PROFILE=1

work_dir="${HOME}/report2_hllc_phase_profile_20260714"
build_dir="${work_dir}/build-report2-hllc-cuda"
export TMPDIR="${work_dir}/.compiler_tmp"
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

echo "SMOKE_TMPDIR=${TMPDIR}"
df -h "${TMPDIR}" | tail -n 1

cmake -S "${work_dir}" -B "${build_dir}" \
    -DFP_ENABLE_CUDA=ON \
    -DFP_ENABLE_OPENMP=ON \
    -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.9/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "${build_dir}" -j 8

binary="${build_dir}/final_project"
if [[ ! -x "${binary}" ]]; then
    echo "PHASE_PROFILE_SMOKE_FAIL: missing executable ${binary}" >&2
    exit 2
fi

set +e
output="$(${binary} \
    --backend=cuda \
    --problem=gaussian_advection_2d \
    --method=weno2 \
    --nx=25 \
    --ny=25 \
    --boundary=periodic \
    --final_time=0.002 \
    --cfl=0.4 \
    --output=none 2>&1)"
solver_status=$?
set -e

printf '%s\n' "${output}"
echo "PHASE_PROFILE_SOLVER_EXIT=${solver_status}"
if [[ ${solver_status} -ne 0 ]]; then
    echo "PHASE_PROFILE_SMOKE_FAIL: solver returned ${solver_status}" >&2
    exit "${solver_status}"
fi
grep -q 'CUDA compiled. Detected device count: 1' <<<"${output}"
grep -q 'cuda smoke test : passed' <<<"${output}"
grep -q '^PHASE_PROFILE ' <<<"${output}"
grep -q 'cuda 2D weno2 solver : passed' <<<"${output}"

echo "PHASE_PROFILE_SMOKE_PASS"
'@

$remote = $remote.Replace([string][char]13, "")
$remote | ssh "${User}@${HostName}" bash -s
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC phase-profile smoke test failed with exit code $LASTEXITCODE"
}
