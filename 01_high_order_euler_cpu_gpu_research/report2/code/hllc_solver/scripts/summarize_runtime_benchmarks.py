#!/usr/bin/env python3
"""Print a compact terminal summary of runtime_benchmarks.csv."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    with args.csv_path.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit("benchmark csv is empty")

    serial_best: dict[tuple[str, str], float] = {}
    for row in rows:
        if row["backend"] == "serial":
            serial_best[(row["case"], row["method"])] = float(row["best_seconds"])

    current_case = None
    for row in rows:
        case = row["case"]
        if case != current_case:
            current_case = case
            print(f"\ncase: {case} ({row['problem']})")
            print("backend\tmethod\tbest(s)\tmean(s)\tmedian(s)\tspeedup\tnx\tny")

        best = float(row["best_seconds"])
        baseline = serial_best.get((row["case"], row["method"]))
        speedup = baseline / best if baseline and best > 0.0 else 0.0

        print(
            f"{row['backend']}\t"
            f"{row['method']}\t"
            f"{best:.6f}\t"
            f"{float(row['mean_seconds']):.6f}\t"
            f"{float(row['median_seconds']):.6f}\t"
            f"{speedup:.3f}\t"
            f"{row['nx']}\t"
            f"{row['ny']}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
