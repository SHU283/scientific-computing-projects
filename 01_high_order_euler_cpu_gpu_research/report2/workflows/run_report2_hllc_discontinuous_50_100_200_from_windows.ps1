param(
    [Parameter(Mandatory = $true)]
    [string]$User,
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"

$Root = "D:\code\Final_projects"
$Source = Join-Path $Root "report2-hllc"
$Archive = Join-Path $Root "report2_hllc_discontinuous_50_100_200_source.tar.gz"
$RemoteScript = Join-Path $PSScriptRoot "run_report2_hllc_discontinuous_50_100_200_remote.sh"

if (Test-Path -LiteralPath $Archive) {
    Remove-Item -LiteralPath $Archive -Force
}

Write-Host "=== Packing Report 2 HLLC discontinuous 50/100/200 source ==="
tar.exe -czf $Archive -C $Source CMakeLists.txt include src scripts
if ($LASTEXITCODE -ne 0) {
    throw "Source packing failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Uploading source and HLLC discontinuous runner ==="
scp $Archive "${User}@${HostName}:~/report2_hllc_discontinuous_50_100_200_source.tar.gz"
if ($LASTEXITCODE -ne 0) {
    throw "Source upload failed with exit code $LASTEXITCODE"
}
scp $RemoteScript "${User}@${HostName}:~/run_report2_hllc_discontinuous_50_100_200_remote.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Runner upload failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "=== Starting background HLLC discontinuous timing at n=50,100,200 ==="
ssh "${User}@${HostName}" "sed -i 's/\r$//' ~/run_report2_hllc_discontinuous_50_100_200_remote.sh && rm -f ~/report2_hllc_discontinuous_50_100_200.log && nohup bash ~/run_report2_hllc_discontinuous_50_100_200_remote.sh > ~/report2_hllc_discontinuous_50_100_200.log 2>&1 < /dev/null & echo `$!"
if ($LASTEXITCODE -ne 0) {
    throw "Remote HLLC discontinuous launch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Started. Check progress with:"
Write-Host "ssh $User@$HostName `"tail -n 80 ~/report2_hllc_discontinuous_50_100_200.log`""
Write-Host ""
Write-Host "Fetch later with:"
Write-Host ".\fetch_report2_hllc_discontinuous_50_100_200_results.ps1 -User $User -HostName $HostName"
