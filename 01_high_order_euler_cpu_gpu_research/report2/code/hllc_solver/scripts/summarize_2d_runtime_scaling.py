#!/usr/bin/env python3
"""Summarise runtime_scaling_2d.csv by grid resolution."""

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
        raise SystemExit("runtime scaling csv is empty")

    serial_best: dict[tuple[str, int], float] = {}
    for row in rows:
        if row["backend"] == "serial":
            serial_best[(row["method"], int(row["nx"]))] = float(row["best_seconds"])

    current_method = None
    for row in rows:
        method = row["method"]
        if method != current_method:
            current_method = method
            print(f"\nmethod: {method}")
            print("backend\tnx\tbest(s)\tmean(s)\tmedian(s)\tspeedup")

        best = float(row["best_seconds"])
        baseline = serial_best.get((method, int(row["nx"])))
        speedup = baseline / best if baseline and best > 0.0 else 0.0
        print(
            f"{row['backend']}\t"
            f"{row['nx']}\t"
            f"{best:.6f}\t"
            f"{float(row['mean_seconds']):.6f}\t"
            f"{float(row['median_seconds']):.6f}\t"
            f"{speedup:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
