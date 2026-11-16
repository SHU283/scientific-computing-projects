param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Project = Join-Path $Root "report2-hllc"
$RemoteArchive = "~/report2_hllc_1d_results.tar.gz"
$LocalArchive = Join-Path $Root "report2_hllc_1d_results.tar.gz"

Write-Host "=== Checking remote HLLC 1D archive ==="
ssh "${User}@${HostName}" "ls -lh ~/report2_hllc_1d_results.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Downloading HLLC 1D archive ==="
scp "${User}@${HostName}:$RemoteArchive" $LocalArchive
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Updating canonical Report 2 HLLC 1D results ==="
tar.exe -xzf $LocalArchive -C $Project
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Push-Location $Project
try {
    python scripts/plot_report2_convergence_1d.py `
        results/report2_hllc/validation_1d/report2_high_order_convergence_1d.csv
    python scripts/plot_report2_high_resolution_fields.py `
        --base-dir results/report2_hllc/high_resolution_fields `
        --resolution-1d 800 `
        --resolution-2d 600
} finally {
    Pop-Location
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $LocalArchive
