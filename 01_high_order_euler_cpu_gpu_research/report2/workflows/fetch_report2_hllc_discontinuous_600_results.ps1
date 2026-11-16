param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$ResultArchive = Join-Path $Root "report2_hllc_discontinuous_600_results.tar.gz"
$RemoteArchive = "~/report2_hllc_discontinuous_600_results.tar.gz"

Write-Host "=== Checking remote HLLC discontinuous 600 archive ==="
ssh "${User}@${HostName}" "ls -lh ${RemoteArchive}"
if ($LASTEXITCODE -ne 0) {
    throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Fetching HLLC discontinuous 600 results ==="
scp "${User}@${HostName}:${RemoteArchive}" $ResultArchive
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
