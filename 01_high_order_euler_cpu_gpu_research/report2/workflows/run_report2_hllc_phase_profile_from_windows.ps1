param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_phase_profile_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_phase_profile_remote.sh"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing Report 2 HLLC phase-profiling source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts docs/phase_profiling_2d.md
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and phase-profiling runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_phase_profile_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_phase_profile_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting background HLLC phase-profiling campaign ==="
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_phase_profile_remote.sh; rm -f ~/report2_hllc_phase_profile.log; nohup bash ~/run_report2_hllc_phase_profile_remote.sh > ~/report2_hllc_phase_profile.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote background launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Background run started. Check progress with:"
Write-Host "$PSScriptRoot\check_report2_hllc_phase_profile_status.ps1 -User $User -HostName $HostName"
Write-Host ""
Write-Host "Fetch later with:"
Write-Host "$PSScriptRoot\fetch_report2_hllc_phase_profile_results.ps1 -User $User -HostName $HostName"
