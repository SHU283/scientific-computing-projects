param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$archive = Join-Path $root "report2_hllc_gpu_limit_extension_source.tar.gz"
$remoteScript = Join-Path $PSScriptRoot "run_report2_hllc_gpu_limit_extension_remote.sh"

if (Test-Path -LiteralPath $archive) {
    Remove-Item -LiteralPath $archive -Force
}

Write-Host "=== Packing HLLC CUDA-only GPU-limit extension source ==="
$resultDir = Join-Path $source "results\report2_hllc\performance\gpu_limit"
if (Test-Path -LiteralPath $resultDir) {
    tar.exe -czf $archive -C $source CMakeLists.txt include src scripts results/report2_hllc/performance/gpu_limit
} else {
    tar.exe -czf $archive -C $source CMakeLists.txt include src scripts
}
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and HLLC GPU-limit extension runner ==="
scp $archive "${User}@${HostName}:~/report2_hllc_gpu_limit_extension_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $remoteScript "${User}@${HostName}:~/run_report2_hllc_gpu_limit_extension_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting HLLC CUDA-only 2400/3000/3600/4200/4600/5000 extension in the remote background ==="
Write-Host "This is an exploratory one-repeat sweep; each completed case is written to CSV immediately."
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_gpu_limit_extension_remote.sh && rm -f ~/report2_hllc_gpu_limit_extension.log && nohup bash ~/run_report2_hllc_gpu_limit_extension_remote.sh > ~/report2_hllc_gpu_limit_extension.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC GPU-limit extension background launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Started. Check progress with:"
Write-Host "ssh ${User}@${HostName} `"tail -n 80 ~/report2_hllc_gpu_limit_extension.log`""
Write-Host ""
Write-Host "Fetch partial or final results with:"
Write-Host ".\fetch_report2_hllc_gpu_limit_extension_results.ps1 -User ${User} -HostName ${HostName}"
