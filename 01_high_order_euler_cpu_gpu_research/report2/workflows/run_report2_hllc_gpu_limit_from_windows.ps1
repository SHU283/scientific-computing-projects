param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$archive = Join-Path $root "report2_hllc_gpu_limit_source.tar.gz"
$remoteScript = Join-Path $PSScriptRoot "run_report2_hllc_gpu_limit_remote.sh"

if (Test-Path -LiteralPath $archive) {
    Remove-Item -LiteralPath $archive -Force
}

Write-Host "=== Packing HLLC CUDA-only GPU-limit source ==="
tar.exe -czf $archive -C $source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and HLLC GPU-limit runner ==="
scp $archive "${User}@${HostName}:~/report2_hllc_gpu_limit_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $remoteScript "${User}@${HostName}:~/run_report2_hllc_gpu_limit_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting HLLC CUDA-only 600/1000/1400/1800 smooth sweep in the remote background ==="
Write-Host "This tests only the CUDA backend for the official HLLC Report 2 route."
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_gpu_limit_remote.sh && rm -f ~/report2_hllc_gpu_limit.log && nohup bash ~/run_report2_hllc_gpu_limit_remote.sh > ~/report2_hllc_gpu_limit.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC GPU-limit background launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Started. Check progress with:"
Write-Host "ssh ${User}@${HostName} `"tail -n 60 ~/report2_hllc_gpu_limit.log`""
Write-Host ""
Write-Host "Fetch later with:"
Write-Host ".\fetch_report2_hllc_gpu_limit_results.ps1 -User ${User} -HostName ${HostName}"
