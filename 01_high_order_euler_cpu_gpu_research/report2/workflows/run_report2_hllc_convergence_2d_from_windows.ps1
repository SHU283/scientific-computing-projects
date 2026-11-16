param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$archive = Join-Path $root "report2_hllc_convergence_2d_source.tar.gz"
$remoteScript = Join-Path $PSScriptRoot "run_report2_hllc_convergence_2d_remote.sh"
$resultArchive = Join-Path $root "report2_hllc_convergence_2d_results.tar.gz"

if (Test-Path -LiteralPath $archive) {
    Remove-Item -LiteralPath $archive -Force
}

Write-Host "=== Packing Report 2 HLLC 2D convergence source ==="
tar.exe -czf $archive -C $source CMakeLists.txt include src scripts

Write-Host ""
Write-Host "=== Uploading source and HLLC convergence runner ==="
scp $archive "${User}@${HostName}:~/report2_hllc_convergence_2d_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $remoteScript "${User}@${HostName}:~/run_report2_hllc_convergence_2d_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Running HLLC 2D convergence on Lovelace ==="
ssh "${User}@${HostName}" `
    "sed -i 's/\r$//' ~/run_report2_hllc_convergence_2d_remote.sh && bash ~/run_report2_hllc_convergence_2d_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC 2D convergence run failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Fetching HLLC 2D convergence results ==="
scp "${User}@${HostName}:~/report2_hllc_convergence_2d_results.tar.gz" $resultArchive
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Updating canonical HLLC results ==="
tar.exe -xzf $resultArchive -C $source
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $resultArchive
