#!/usr/bin/env python3
"""
Read a 2D discontinuous runtime CSV and print a compact summary table.
"""

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
        raise SystemExit("runtime csv is empty")

    serial_best: dict[tuple[str, str], float] = {}
    for row in rows:
        if row["backend"] == "serial":
            serial_best[(row["case"], row["method"])] = float(row["best_seconds"])

    current_case = None
    current_method = None
    for row in rows:
        case = row["case"]
        method = row["method"]
        if case != current_case or method != current_method:
            current_case = case
            current_method = method
            print(f"\ncase: {case} | method: {method}")
            print("backend\tbest(s)\tmean(s)\tmedian(s)\tspeedup")

        best = float(row["best_seconds"])
        baseline = serial_best.get((case, method))
        speedup = baseline / best if baseline and best > 0.0 else 0.0
        print(
            f"{row['backend']}\t"
            f"{best:.6f}\t"
            f"{float(row['mean_seconds']):.6f}\t"
            f"{float(row['median_seconds']):.6f}\t"
            f"{speedup:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
