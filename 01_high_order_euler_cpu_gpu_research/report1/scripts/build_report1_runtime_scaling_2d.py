#!/usr/bin/env python3
"""
Merge baseline 2D runtime scaling rows with WENO-2 rows for Report 1.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", default="results/stageC/runtime_scaling_2d.csv")
    parser.add_argument("--weno2", default="results/stageC/runtime_scaling_2d_weno2.csv")
    parser.add_argument("--output", default="results/stageC/runtime_scaling_2d_report1.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    rows = load_rows(root / args.baseline) + load_rows(root / args.weno2)

    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "problem",
                "backend",
                "method",
                "nx",
                "ny",
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved merged runtime scaling csv: {output_path}")
    print(f"rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

