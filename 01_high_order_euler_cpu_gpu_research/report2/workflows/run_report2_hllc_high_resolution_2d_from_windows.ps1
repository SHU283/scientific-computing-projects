param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74",
    [int]$Resolution = 600,
    [switch]$Detached
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_high_resolution_2d_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_high_resolution_2d_remote.sh"
$ResultArchive = Join-Path $Root "report2_hllc_high_resolution_2d_results.tar.gz"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing Report 2 HLLC high-resolution 2D source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and HLLC high-resolution runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_high_resolution_2d_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_high_resolution_2d_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Running $Resolution x $Resolution HLLC 2D field exports on Lovelace ==="
if ($Detached) {
    ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_high_resolution_2d_remote.sh && rm -f ~/report2_hllc_high_resolution_2d.log && nohup env REPORT2_HLLC_FIELD_RESOLUTION=$Resolution bash ~/run_report2_hllc_high_resolution_2d_remote.sh > ~/report2_hllc_high_resolution_2d.log 2>&1 < /dev/null & echo `$!"
} else {
    ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_high_resolution_2d_remote.sh && REPORT2_HLLC_FIELD_RESOLUTION=$Resolution bash ~/run_report2_hllc_high_resolution_2d_remote.sh"
}
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC high-resolution field run failed with exit code $LASTEXITCODE"
}

if ($Detached) {
    Write-Host ""
    Write-Host "Detached Lovelace job started."
    Write-Host "Remote log: ~/report2_hllc_high_resolution_2d.log"
    Write-Host "Fetch later with:"
    Write-Host ".\fetch_report2_hllc_high_resolution_2d_results.ps1 -User $User -HostName $HostName"
    exit 0
}

Write-Host ""
Write-Host "=== Fetching HLLC high-resolution field results ==="
scp "${User}@${HostName}:~/report2_hllc_high_resolution_2d_results.tar.gz" $ResultArchive
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
