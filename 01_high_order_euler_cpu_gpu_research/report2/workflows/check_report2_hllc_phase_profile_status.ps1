param(
    [string]$User = "sx284",
    [string]$HostName = "172.24.62.74"
)

$ErrorActionPreference = "Stop"
$remote = @'
echo STATUS_BEGIN
date
echo "=== active processes ==="
pgrep -af '[r]eport2_hllc_phase_profile' || true
pgrep -af '[f]inal_project.*report2_hllc_phase_profile' || true
echo "=== archive ==="
ls -lh ~/report2_hllc_phase_profile_results.tar.gz 2>/dev/null || echo NO_ARCHIVE_YET
echo "=== CSV ==="
csv=~/report2_hllc_phase_profile_20260714/results/report2_hllc/performance/phase_profile/phase_profile_2d_hllc_lovelace_a30.csv
if [ -f "$csv" ]; then
    wc -l "$csv"
    tail -n 8 "$csv"
else
    echo NO_CSV_YET
fi
echo "=== log tail ==="
tail -n 40 ~/report2_hllc_phase_profile.log 2>/dev/null || true
echo STATUS_END
'@
$remote = $remote.Replace([string][char]13, "")
$remote | ssh "${User}@${HostName}" bash -s
if ($LASTEXITCODE -ne 0) {
    throw "Remote status check failed with exit code $LASTEXITCODE"
}
