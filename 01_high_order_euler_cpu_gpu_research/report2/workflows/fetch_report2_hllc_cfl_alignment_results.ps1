param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$root = "D:\code\Final_projects"
$project = Join-Path $root "report2-hllc"
$archiveName = "report2_hllc_cfl_alignment_results.tar.gz"
$localArchive = Join-Path $root $archiveName

Write-Host "=== Fetching CFL-aligned HLLC validation results ==="
scp "${User}@${HostName}:~/$archiveName" $localArchive
if ($LASTEXITCODE -ne 0) { throw "Result download failed with exit code $LASTEXITCODE" }

tar.exe -xzf $localArchive -C $project
if ($LASTEXITCODE -ne 0) { throw "Result extraction failed with exit code $LASTEXITCODE" }

python (Join-Path $project "scripts\audit_report2_cfl_policy.py")
if ($LASTEXITCODE -ne 0) { throw "CFL audit failed with exit code $LASTEXITCODE" }

Write-Host "Done: $localArchive"
