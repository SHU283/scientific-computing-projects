param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_openmp_thread_scaling_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_openmp_thread_scaling_remote.sh"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing Report 2 HLLC OpenMP thread-scaling source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and remote runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_openmp_thread_scaling_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_openmp_thread_scaling_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting background HLLC OpenMP thread-scaling run ==="
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_openmp_thread_scaling_remote.sh; rm -f ~/report2_hllc_openmp_thread_scaling.log; nohup bash ~/run_report2_hllc_openmp_thread_scaling_remote.sh > ~/report2_hllc_openmp_thread_scaling.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote background launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Background run started. Check progress with:"
Write-Host "ssh ${User}@${HostName} `"tail -n 60 ~/report2_hllc_openmp_thread_scaling.log`""
Write-Host ""
Write-Host "Fetch later with:"
Write-Host "$PSScriptRoot\fetch_report2_hllc_openmp_thread_scaling_results.ps1 -User $User -HostName $HostName"
