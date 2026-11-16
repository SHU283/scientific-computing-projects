#!/usr/bin/env python3
"""Audit completeness and timing consistency of Report 2 phase profiles."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
BACKENDS = ("serial", "openmp", "cuda")
EXPECTED_REPEATS = {100: 3, 600: 1}
RAW_PHASES = (
    "setup_s",
    "cfl_s",
    "workspace_s",
    "flux_x_s",
    "flux_y_s",
    "spatial_fused_s",
    "rhs_assembly_s",
    "state_update_s",
    "admissibility_s",
    "transfer_s",
    "unattributed_s",
)


def as_float(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {row}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    with args.csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("phase-profile CSV is empty")

    output_dir = args.output_dir or args.csv_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    issues: list[str] = []

    for row_number, row in enumerate(rows, start=2):
        key = (row["method"], row["backend"], int(row["nx"]))
        grouped[key].append(row)
        if int(row["nx"]) != int(row["ny"]):
            issues.append(f"row {row_number}: non-square grid")
        if row["problem"] != "gaussian_advection_2d":
            issues.append(f"row {row_number}: unexpected problem {row['problem']}")
        if int(row["retries"]) != 0:
            issues.append(f"row {row_number}: retries={row['retries']}")
        total = as_float(row, "total_s")
        wall = as_float(row, "wall_seconds")
        phases = sum(as_float(row, field) for field in RAW_PHASES)
        closure = abs(phases - total) / max(total, 1.0e-15)
        if closure > 5.0e-3:
            issues.append(f"row {row_number}: phase closure error {closure:.2%}")
        # wall_seconds is measured by the Python launcher and includes process
        # startup plus CUDA context initialisation.  It is intentionally not
        # expected to close against the solver-internal phase timer.
        if any(as_float(row, field) < 0.0 for field in ("total_s", *RAW_PHASES)):
            issues.append(f"row {row_number}: negative timing")

    expected_keys = {
        (method, backend, resolution)
        for resolution in EXPECTED_REPEATS
        for method in METHODS
        for backend in BACKENDS
    }
    for key in sorted(expected_keys - set(grouped)):
        issues.append(f"missing group: method={key[0]} backend={key[1]} n={key[2]}")
    for key, group in sorted(grouped.items()):
        expected = EXPECTED_REPEATS.get(key[2])
        if expected is None:
            issues.append(f"unexpected resolution in group {key}")
        elif len(group) != expected:
            issues.append(f"group {key}: expected {expected} rows, found {len(group)}")

    audit_path = output_dir / "phase_profile_2d_hllc.audit.csv"
    fields = (
        "method",
        "backend",
        "resolution",
        "repeats",
        "median_total_s",
        "min_total_s",
        "max_total_s",
        "relative_range",
        "median_wall_over_total",
        "median_process_overhead_s",
        "steps_min",
        "steps_max",
        "retries_sum",
        "median_accounted_fraction",
    )
    audit_rows: list[dict[str, object]] = []
    for (method, backend, resolution), group in sorted(
        grouped.items(), key=lambda item: (item[0][2], METHODS.index(item[0][0]), BACKENDS.index(item[0][1]))
    ):
        totals = [as_float(row, "total_s") for row in group]
        walls = [as_float(row, "wall_seconds") for row in group]
        accounted = [as_float(row, "accounted_s") for row in group]
        median_total = statistics.median(totals)
        audit_rows.append(
            {
                "method": method,
                "backend": backend,
                "resolution": resolution,
                "repeats": len(group),
                "median_total_s": median_total,
                "min_total_s": min(totals),
                "max_total_s": max(totals),
                "relative_range": (max(totals) - min(totals)) / max(median_total, 1.0e-15),
                "median_wall_over_total": statistics.median(
                    wall / max(total, 1.0e-15) for wall, total in zip(walls, totals)
                ),
                "median_process_overhead_s": statistics.median(
                    max(0.0, wall - total) for wall, total in zip(walls, totals)
                ),
                "steps_min": min(int(row["steps"]) for row in group),
                "steps_max": max(int(row["steps"]) for row in group),
                "retries_sum": sum(int(row["retries"]) for row in group),
                "median_accounted_fraction": statistics.median(
                    value / max(total, 1.0e-15) for value, total in zip(accounted, totals)
                ),
            }
        )

    with audit_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(audit_rows)

    report_path = output_dir / "phase_profile_2d_hllc.audit.md"
    lines = [
        "# HLLC phase-profile audit",
        "",
        f"- Raw rows: {len(rows)}",
        f"- Expected rows: {sum(EXPECTED_REPEATS.values()) * len(METHODS) * len(BACKENDS)}",
        f"- Audit findings: {len(issues)}",
        "",
        "| n | method | backend | repeats | median total (s) | range / median | accounted | process overhead (s) |",
        "|---:|:---|:---|---:|---:|---:|---:|---:|",
    ]
    for row in audit_rows:
        lines.append(
            f"| {row['resolution']} | {str(row['method']).upper()} | {row['backend']} | "
            f"{row['repeats']} | {float(row['median_total_s']):.6g} | "
            f"{100.0 * float(row['relative_range']):.2f}% | "
            f"{100.0 * float(row['median_accounted_fraction']):.2f}% | "
            f"{float(row['median_process_overhead_s']):.4g} |"
        )
    lines.extend(["", "## Findings", ""])
    lines.extend(f"- {issue}" for issue in issues)
    if not issues:
        lines.append("- No completeness or timing-consistency problems detected.")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"saved audit: {audit_path}")
    print(f"saved audit report: {report_path}")
    print(f"audit findings: {len(issues)}")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
