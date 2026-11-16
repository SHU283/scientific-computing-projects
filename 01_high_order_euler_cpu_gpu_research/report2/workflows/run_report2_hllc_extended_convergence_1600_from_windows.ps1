param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$archiveName = "report2_hllc_extended_convergence_1600_source.tar.gz"
$archive = Join-Path $root $archiveName
$remoteScript = Join-Path $PSScriptRoot "run_report2_hllc_extended_convergence_1600_remote.sh"

if (Test-Path -LiteralPath $archive) {
    Remove-Item -LiteralPath $archive -Force
}

Write-Host "=== Packing HLLC n=1600 convergence extension ==="
tar.exe -czf $archive -C $source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) { throw "Source packing failed with exit code $LASTEXITCODE" }

Write-Host "=== Uploading source and remote runner ==="
scp $archive "${User}@${HostName}:~/$archiveName"
if ($LASTEXITCODE -ne 0) { throw "Source upload failed with exit code $LASTEXITCODE" }
scp $remoteScript "${User}@${HostName}:~/run_report2_hllc_extended_convergence_1600_remote.sh"
if ($LASTEXITCODE -ne 0) { throw "Runner upload failed with exit code $LASTEXITCODE" }

Write-Host "=== Starting background n=1600 convergence extension ==="
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_extended_convergence_1600_remote.sh && rm -f ~/report2_hllc_extended_convergence_1600.log && nohup bash ~/run_report2_hllc_extended_convergence_1600_remote.sh > ~/report2_hllc_extended_convergence_1600.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) { throw "Remote launch failed with exit code $LASTEXITCODE" }

Write-Host "Background run started. Check with:"
Write-Host "  $PSScriptRoot\check_report2_hllc_extended_convergence_1600_status.ps1 -User $User -HostName $HostName"
Write-Host "Fetch later with:"
Write-Host "  $PSScriptRoot\fetch_report2_hllc_extended_convergence_1600_results.ps1 -User $User -HostName $HostName"
