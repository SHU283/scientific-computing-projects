param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$root = "D:\code\Final_projects"
$project = Join-Path $root "report2-hllc"
$archiveName = "report2_hllc_extended_convergence_1600_results.tar.gz"
$localArchive = Join-Path $root $archiveName

Write-Host "=== Checking remote n=1600 convergence archive ==="
ssh "${User}@${HostName}" "ls -lh ~/$archiveName"
if ($LASTEXITCODE -ne 0) { throw "Remote archive was not found or SSH failed with exit code $LASTEXITCODE" }

Write-Host "=== Downloading n=1600 convergence results ==="
scp "${User}@${HostName}:~/$archiveName" $localArchive
if ($LASTEXITCODE -ne 0) { throw "Result download failed with exit code $LASTEXITCODE" }

tar.exe -xzf $localArchive -C $project
if ($LASTEXITCODE -ne 0) { throw "Result extraction failed with exit code $LASTEXITCODE" }

Write-Host "Done: $localArchive"
