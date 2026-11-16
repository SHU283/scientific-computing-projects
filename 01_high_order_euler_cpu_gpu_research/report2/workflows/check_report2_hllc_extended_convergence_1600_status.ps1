param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$remoteCommand = @'
echo STATUS_BEGIN
date
echo === active ===
pgrep -af '[r]un_report2_hllc_extended_convergence_1600_remote.sh'
pgrep -af '[f]inal_project.*report2_hllc_extended_convergence_1600'
echo === archive ===
if [ -f "$HOME/report2_hllc_extended_convergence_1600_results.tar.gz" ]; then
    ls -lh "$HOME/report2_hllc_extended_convergence_1600_results.tar.gz"
else
    echo NO_ARCHIVE_YET
fi
echo === partial results ===
find "$HOME/report2_hllc_extended_convergence_1600_20260721/results/report2_hllc/extended_convergence_1600" \
    -type f -maxdepth 3 -printf '%f %s bytes\n' 2>/dev/null | sort
echo === log tail ===
tail -n 50 "$HOME/report2_hllc_extended_convergence_1600.log" 2>/dev/null
echo STATUS_END
'@

ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) { throw "Remote status check failed with exit code $LASTEXITCODE" }
