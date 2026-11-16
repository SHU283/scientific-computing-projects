param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$Root = "D:\code\Final_projects"
$RemoteScript = Join-Path $PSScriptRoot "resume_report2_hllc_phase_profile_remote.sh"
$SmokeScript = Join-Path $Root "smoke_report2_hllc_phase_profile_from_windows.ps1"

Write-Host "=== Running mandatory HLLC phase-profile smoke test ==="
& $SmokeScript -User $User -HostName $HostName
if ($LASTEXITCODE -ne 0) {
    throw "Phase-profile smoke test did not pass; the long resume was not launched."
}

Write-Host "=== Uploading HLLC phase-profile resume runner ==="
scp $RemoteScript "${User}@${HostName}:~/resume_report2_hllc_phase_profile_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Resume runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Resuming HLLC phase profiling with CUDA 12.9 ==="
$remoteCommand = "sed -i 's/\r$//' ~/resume_report2_hllc_phase_profile_remote.sh; rm -f ~/report2_hllc_phase_profile.log; nohup bash ~/resume_report2_hllc_phase_profile_remote.sh > ~/report2_hllc_phase_profile.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote resume launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Background resume started. Check with:"
Write-Host "$PSScriptRoot\check_report2_hllc_phase_profile_status.ps1 -User $User -HostName $HostName"
Write-Host ""
Write-Host "Fetch later with:"
Write-Host "$PSScriptRoot\fetch_report2_hllc_phase_profile_results.ps1 -User $User -HostName $HostName"
