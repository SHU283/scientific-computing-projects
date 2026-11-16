param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$destination = Join-Path $root "report2_hllc_gpu_limit_results.tar.gz"
$remoteArchive = "~/report2_hllc_gpu_limit_results.tar.gz"

Write-Host "=== Checking remote HLLC CUDA-only GPU-limit archive ==="
ssh "${User}@${HostName}" "ls -lh ${remoteArchive}"
if ($LASTEXITCODE -ne 0) {
    throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Downloading HLLC CUDA-only GPU-limit results ==="
scp "${User}@${HostName}:${remoteArchive}" $destination
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

if (-not (Test-Path -LiteralPath $destination)) {
    throw "Result archive was not created locally: $destination"
}

Write-Host ""
Write-Host "=== Updating canonical HLLC GPU-limit results ==="
tar.exe -xzf $destination -C $source
if ($LASTEXITCODE -ne 0) {
    throw "Result extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $destination
