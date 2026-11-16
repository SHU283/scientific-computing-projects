param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_backend_consistency_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_backend_consistency_remote.sh"
$ResultArchive = Join-Path $Root "report2_hllc_backend_consistency_results.tar.gz"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing Report 2 HLLC backend-consistency source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and HLLC backend-consistency runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_backend_consistency_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_backend_consistency_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Running representative Report 2 HLLC backend-consistency checks ==="
ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_backend_consistency_remote.sh && bash ~/run_report2_hllc_backend_consistency_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC backend-consistency campaign failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Fetching HLLC backend-consistency results ==="
scp "${User}@${HostName}:~/report2_hllc_backend_consistency_results.tar.gz" $ResultArchive
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Updating canonical report2-hllc results ==="
tar.exe -xzf $ResultArchive -C $Source
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $ResultArchive
