param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$LocalArchive = Join-Path $Root "report2_hllc_weno2_openmp_refresh_results.tar.gz"
$RemoteArchive = "/home/raid/$User/report2_hllc_weno2_openmp_refresh_results.tar.gz"

Write-Host "=== Checking remote WENO-2 OpenMP refresh archive ==="
ssh "${User}@${HostName}" "ls -lh '$RemoteArchive'"
if ($LASTEXITCODE -ne 0) {
    throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Downloading refresh candidates ==="
scp "${User}@${HostName}:$RemoteArchive" $LocalArchive
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Extracting candidates without replacing canonical CSVs ==="
tar.exe -xzf $LocalArchive -C $Source
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Candidate results are under:"
Write-Host (Join-Path $Source "results\report2_hllc\performance\weno2_openmp_refresh")
