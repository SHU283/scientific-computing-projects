param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_weno2_openmp_refresh_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_weno2_openmp_refresh_remote.sh"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing audited HLLC WENO-2 OpenMP refresh source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and refresh runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_weno2_openmp_refresh_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_weno2_openmp_refresh_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting background WENO-2 OpenMP refresh ==="
$RemoteCommand = "sed -i 's/`r$//' ~/run_report2_hllc_weno2_openmp_refresh_remote.sh; rm -f ~/report2_hllc_weno2_openmp_refresh.log; nohup bash ~/run_report2_hllc_weno2_openmp_refresh_remote.sh > ~/report2_hllc_weno2_openmp_refresh.log 2>&1 < /dev/null & echo `$!"
ssh "${User}@${HostName}" $RemoteCommand
if ($LASTEXITCODE -ne 0) {
    throw "Remote background launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Background refresh started. Check with:"
Write-Host "ssh ${User}@${HostName} `"tail -n 80 ~/report2_hllc_weno2_openmp_refresh.log`""
Write-Host ""
Write-Host "Fetch later with:"
Write-Host "$PSScriptRoot\fetch_report2_hllc_weno2_openmp_refresh_results.ps1 -User $User -HostName $HostName"
