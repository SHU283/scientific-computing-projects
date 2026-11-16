param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Project = Join-Path $Root "report2-hllc"
$RemoteArchive = "~/report2_hllc_openmp_thread_scaling_results.tar.gz"
$LocalArchive = Join-Path $Root "report2_hllc_openmp_thread_scaling_results.tar.gz"

Write-Host "=== Checking remote HLLC OpenMP thread-scaling archive ==="
ssh "${User}@${HostName}" "ls -lh ~/report2_hllc_openmp_thread_scaling_results.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Downloading HLLC OpenMP thread-scaling archive ==="
scp "${User}@${HostName}:$RemoteArchive" $LocalArchive
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}
if (-not (Test-Path -LiteralPath $LocalArchive)) {
    throw "Result download finished but local archive was not created: $LocalArchive"
}
if ((Get-Item -LiteralPath $LocalArchive).Length -le 0) {
    throw "Result download created an empty archive: $LocalArchive"
}

Write-Host ""
Write-Host "=== Updating canonical Report 2 HLLC OpenMP thread-scaling results ==="
tar.exe -xzf $LocalArchive -C $Project
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $LocalArchive
