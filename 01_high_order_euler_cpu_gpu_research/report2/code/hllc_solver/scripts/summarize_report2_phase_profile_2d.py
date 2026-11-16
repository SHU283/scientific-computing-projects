#!/usr/bin/env python3
"""Build compact CSV and Markdown summaries from the HLLC phase profiles."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


BACKENDS = ("serial", "openmp", "cuda")
PHASES = (
    "setup_transfer_s",
    "cfl_s",
    "spatial_operator_s",
    "update_validity_s",
    "other_s",
)


def median_rows(path: Path) -> dict[tuple[str, str, int], dict[str, float]]:
    grouped: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            grouped[(row["method"], row["backend"], int(row["nx"]))].append(row)

    result: dict[tuple[str, str, int], dict[str, float]] = {}
    for key, rows in grouped.items():
        def med(field: str) -> float:
            return statistics.median(float(row[field]) for row in rows)

        result[key] = {
            "total_s": med("total_s"),
            "setup_transfer_s": med("setup_s") + med("workspace_s") + med("transfer_s"),
            "cfl_s": med("cfl_s"),
            "spatial_operator_s": med("spatial_operator_s"),
            "update_validity_s": med("state_update_s") + med("admissibility_s"),
            "other_s": med("unattributed_s"),
        }
    return result


def repeat_variability(path: Path) -> dict[tuple[str, str, int], tuple[int, float]]:
    grouped: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            grouped[(row["method"], row["backend"], int(row["nx"]))].append(
                float(row["total_s"])
            )

    result: dict[tuple[str, str, int], tuple[int, float]] = {}
    for key, values in grouped.items():
        median = statistics.median(values)
        spread = (max(values) - min(values)) / median if median > 0.0 else 0.0
        result[key] = (len(values), spread)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    rows = median_rows(args.csv_path)
    variability = repeat_variability(args.csv_path)
    output_dir = args.output_dir or args.csv_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = output_dir / "phase_profile_2d_hllc.summary.csv"
    fieldnames = (
        "method", "backend", "resolution", "total_s", "total_speedup_vs_serial",
        "dominant_phase", "dominant_fraction", *PHASES,
    )

    summary_rows: list[dict[str, object]] = []
    for method, backend, resolution in sorted(rows, key=lambda key: (key[2], key[0], BACKENDS.index(key[1]))):
        data = rows[(method, backend, resolution)]
        serial_total = rows[(method, "serial", resolution)]["total_s"]
        dominant = max(PHASES, key=lambda phase: data[phase])
        summary_rows.append({
            "method": method,
            "backend": backend,
            "resolution": resolution,
            "total_s": data["total_s"],
            "total_speedup_vs_serial": serial_total / data["total_s"],
            "dominant_phase": dominant.removesuffix("_s"),
            "dominant_fraction": data[dominant] / data["total_s"],
            **{phase: data[phase] for phase in PHASES},
        })

    with summary_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    markdown = output_dir / "phase_profile_2d_hllc.analysis.md"
    lines = [
        "# HLLC 2D phase-profile summary",
        "",
        "Times are medians where repeats are available. Speedups use the profiled serial run",
        "for the same method and resolution. CUDA phase timings include diagnostic boundary",
        "synchronizations and are therefore not replacements for the canonical no-I/O timings.",
        "",
        "| n | method | backend | total (s) | speedup | dominant phase | share |",
        "|---:|:---|:---|---:|---:|:---|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['resolution']} | {str(row['method']).upper()} | {row['backend']} | "
            f"{float(row['total_s']):.4g} | {float(row['total_speedup_vs_serial']):.3f} | "
            f"{row['dominant_phase']} | {100.0 * float(row['dominant_fraction']):.1f}% |"
        )

    large_resolution = max(key[2] for key in rows)
    lines.extend([
        "",
        f"## Large-grid phase composition ({large_resolution} x {large_resolution})",
        "",
        "`setup/workspace/transfer` combines initialization, temporary flux-array allocation,",
        "and explicit transfer timing. It is useful as a cost fraction, but its CPU and CUDA",
        "implementation boundaries are not identical, so it is excluded from phase-speedup claims.",
        "",
        "| method | backend | setup/workspace/transfer | CFL | spatial | update/validity | other |",
        "|:---|:---|---:|---:|---:|---:|---:|",
    ])
    for method in ("weno2", "weno3", "ader3"):
        for backend in BACKENDS:
            data = rows[(method, backend, large_resolution)]
            total = data["total_s"]
            lines.append(
                f"| {method.upper()} | {backend} | "
                f"{100.0 * data['setup_transfer_s'] / total:.1f}% | "
                f"{100.0 * data['cfl_s'] / total:.1f}% | "
                f"{100.0 * data['spatial_operator_s'] / total:.1f}% | "
                f"{100.0 * data['update_validity_s'] / total:.1f}% | "
                f"{100.0 * data['other_s'] / total:.1f}% |"
            )

    lines.extend([
        "",
        "## Comparable large-grid phase speedups",
        "",
        "Ratios below use the profiled serial phase for the same method. They diagnose where",
        "acceleration is obtained; they are not substitutes for canonical no-I/O speedups.",
        "",
        "| method | backend | CFL | spatial operator | update/validity |",
        "|:---|:---|---:|---:|---:|",
    ])
    for method in ("weno2", "weno3", "ader3"):
        serial = rows[(method, "serial", large_resolution)]
        for backend in ("openmp", "cuda"):
            data = rows[(method, backend, large_resolution)]
            lines.append(
                f"| {method.upper()} | {backend} | "
                f"{serial['cfl_s'] / data['cfl_s']:.2f}x | "
                f"{serial['spatial_operator_s'] / data['spatial_operator_s']:.2f}x | "
                f"{serial['update_validity_s'] / data['update_validity_s']:.2f}x |"
            )

    unstable = sorted(
        (key, count, spread)
        for key, (count, spread) in variability.items()
        if count > 1 and spread >= 0.20
    )
    lines.extend([
        "",
        "## Interpretation and data-quality notes",
        "",
        f"- The spatial operator is the dominant phase for every {large_resolution} x "
        f"{large_resolution} run. Its share remains especially high for ADER-3, so further "
        "ADER-3 optimization should start with predictor/flux kernels rather than peripheral work.",
        "- OpenMP reduces the spatial cost strongly, but WENO-2 and WENO-3 retain substantial "
        "workspace, CFL, and update overhead. This explains why their end-to-end OpenMP gains "
        "are smaller than their spatial-operator gains.",
        "- CUDA also leaves the spatial operator dominant at large resolution. Kernel and memory-"
        "access optimization therefore remain more valuable than optimizing final diagnostics.",
        f"- The {large_resolution} x {large_resolution} rows have one repeat each, so their phase "
        "fractions are diagnostic rather than uncertainty estimates.",
    ])
    if unstable:
        lines.append("- Small-grid groups with at least 20% range/median variation are:")
        for (method, backend, resolution), count, spread in unstable:
            lines.append(
                f"  - {method.upper()} {backend} at {resolution} x {resolution}: "
                f"{count} repeats, {100.0 * spread:.1f}% range/median."
            )
        lines.append(
            "  These cases should illustrate startup/profiling sensitivity, not supply final "
            "performance claims."
        )
    markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"saved summary: {summary_csv}")
    print(f"saved analysis: {markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
