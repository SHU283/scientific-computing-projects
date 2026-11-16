param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$remoteCommand = @'
echo STATUS_BEGIN
date
echo === active ===
pgrep -af '[r]un_report2_hllc_cfl_alignment_remote.sh'
pgrep -af '[f]inal_project.*report2_hllc_cfl_alignment'
echo === archive ===
if [ -f "$HOME/report2_hllc_cfl_alignment_results.tar.gz" ]; then
    ls -lh "$HOME/report2_hllc_cfl_alignment_results.tar.gz"
else
    echo NO_ARCHIVE_YET
fi
echo === log tail ===
tail -n 40 "$HOME/report2_hllc_cfl_alignment.log" 2>/dev/null
echo STATUS_END
'@

ssh "${User}@${HostName}" $remoteCommand
if ($LASTEXITCODE -ne 0) { throw "Remote status check failed with exit code $LASTEXITCODE" }
