param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$root = "D:\code\Final_projects"
$source = Join-Path $root "report2-hllc"
$archive = Join-Path $root "report2_hllc_cfl_alignment_source.tar.gz"
$remoteScript = Join-Path $PSScriptRoot "run_report2_hllc_cfl_alignment_remote.sh"

if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive -Force }

Write-Host "=== Packing CFL-aligned HLLC validation source ==="
tar.exe -czf $archive -C $source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) { throw "Source packing failed with exit code $LASTEXITCODE" }

Write-Host "=== Uploading source and remote runner ==="
scp $archive "${User}@${HostName}:~/report2_hllc_cfl_alignment_source.tar.gz"
if ($LASTEXITCODE -ne 0) { throw "Source upload failed with exit code $LASTEXITCODE" }
scp $remoteScript "${User}@${HostName}:~/run_report2_hllc_cfl_alignment_remote.sh"
if ($LASTEXITCODE -ne 0) { throw "Runner upload failed with exit code $LASTEXITCODE" }

Write-Host "=== Starting consolidated CFL-aligned validation in the background ==="
$remoteCommand = "sed -i 's/\r$//' ~/run_report2_hllc_cfl_alignment_remote.sh && rm -f ~/report2_hllc_cfl_alignment.log && nohup bash ~/run_report2_hllc_cfl_alignment_remote.sh > ~/report2_hllc_cfl_alignment.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) { throw "Remote campaign launch failed with exit code $LASTEXITCODE" }

Write-Host "Background run started. Check later with:"
Write-Host "  $PSScriptRoot\check_report2_hllc_cfl_alignment_status.ps1 -User $User -HostName $HostName"
Write-Host "Fetch when complete with:"
Write-Host "  $PSScriptRoot\fetch_report2_hllc_cfl_alignment_results.ps1 -User $User -HostName $HostName"
