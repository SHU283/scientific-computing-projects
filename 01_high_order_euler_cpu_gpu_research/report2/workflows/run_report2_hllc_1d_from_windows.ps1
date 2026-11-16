param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$LocalPackage = Join-Path $Root "report2_hllc_1d_source.tar.gz"
$RemotePackage = "~/report2_hllc_1d_source.tar.gz"
$RemoteRunner = "~/run_report2_hllc_1d_remote.sh"
$Runner = Join-Path $PSScriptRoot "run_report2_hllc_1d_remote.sh"

Write-Host "=== Packing Report 2 HLLC 1D source ==="
if (Test-Path -LiteralPath $LocalPackage) {
    Remove-Item -LiteralPath $LocalPackage -Force
}

tar.exe -czf $LocalPackage -C $Source CMakeLists.txt include src scripts

Write-Host ""
Write-Host "=== Uploading source and HLLC 1D runner ==="
scp $LocalPackage "${User}@${HostName}:$RemotePackage"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}

scp $Runner "${User}@${HostName}:$RemoteRunner"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Running HLLC 1D convergence and field exports on Lovelace ==="
ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_1d_remote.sh && bash ~/run_report2_hllc_1d_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC 1D campaign failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Fetching HLLC 1D results ==="
& (Join-Path $Root "fetch_report2_hllc_1d_results.ps1") -User $User -HostName $HostName
