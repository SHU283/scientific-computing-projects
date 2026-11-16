param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74",
    [int]$Resolution = 600
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_discontinuous_evolution_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_discontinuous_evolution_remote.sh"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing HLLC discontinuous-evolution source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and evolution runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_discontinuous_evolution_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_discontinuous_evolution_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting background HLLC discontinuous-evolution snapshots ==="
ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_discontinuous_evolution_remote.sh && rm -f ~/report2_hllc_discontinuous_evolution.log && nohup env REPORT2_HLLC_EVOLUTION_RESOLUTION=$Resolution bash ~/run_report2_hllc_discontinuous_evolution_remote.sh > ~/report2_hllc_discontinuous_evolution.log 2>&1 < /dev/null & echo `$!"
if ($LASTEXITCODE -ne 0) {
    throw "Remote evolution launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Background run started. Check with:"
Write-Host "ssh ${User}@${HostName} `"tail -n 80 ~/report2_hllc_discontinuous_evolution.log`""
Write-Host ""
Write-Host "Fetch later with:"
Write-Host ".\fetch_report2_hllc_discontinuous_evolution_results.ps1 -User $User -HostName $HostName"
