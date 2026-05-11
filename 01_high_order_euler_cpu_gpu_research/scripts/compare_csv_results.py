#!/usr/bin/env python3
"""
Compare two solver output CSV files field by field.

Supported layouts:
- 1D: x,rho,u,p,energy
- 2D: i,j,rho,u,v,p,energy
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
    parser.add_argument("csv_a", type=Path)
    parser.add_argument("csv_b", type=Path)
    args = parser.parse_args()

    if not args.csv_a.exists():
        raise SystemExit(f"file not found: {args.csv_a}")
    if not args.csv_b.exists():
        raise SystemExit(f"file not found: {args.csv_b}")

    rows_a = load_rows(args.csv_a)
    rows_b = load_rows(args.csv_b)

    if len(rows_a) != len(rows_b):
        raise SystemExit(f"row count mismatch: {len(rows_a)} vs {len(rows_b)}")

    keys = [k for k in rows_a[0].keys() if k not in {"x", "i", "j"}]
    max_diff = {k: 0.0 for k in keys}

    for row_a, row_b in zip(rows_a, rows_b):
        for key in keys:
            diff = abs(float(row_a[key]) - float(row_b[key]))
            if diff > max_diff[key]:
                max_diff[key] = diff

    print(f"compare: {args.csv_a.name} vs {args.csv_b.name}")
    for key in keys:
        print(f"max |delta {key}| = {max_diff[key]:.12e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
