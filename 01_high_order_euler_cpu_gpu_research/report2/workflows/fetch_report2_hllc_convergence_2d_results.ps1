param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Project = Join-Path $Root "report2-hllc"
$ArchiveName = "report2_hllc_convergence_2d_results.tar.gz"
$LocalArchive = Join-Path $Root $ArchiveName
$RemoteArchive = "~/$ArchiveName"

Write-Host ""
Write-Host "=== Fetching Report 2 HLLC 2D convergence results ==="
Write-Host "Remote: ${User}@${HostName}:$RemoteArchive"
Write-Host "Local : $LocalArchive"
Write-Host "You may be asked for the Lovelace password."

scp "${User}@${HostName}:$RemoteArchive" "$LocalArchive"
if ($LASTEXITCODE -ne 0) {
    throw "Result download failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Updating canonical report2-hllc results ==="
tar -xzf "$LocalArchive" -C "$Project"
if ($LASTEXITCODE -ne 0) {
    throw "Archive extraction failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Done. Result archive:"
Write-Host $LocalArchive
Write-Host ""
Write-Host "Updated result files under:"
Write-Host (Join-Path $Project "results\report2_hllc\validation_2d_lovelace")
