#!/usr/bin/env python3
"""
Print a compact summary from weno3_sensitivity_1d.csv.
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
        raise SystemExit("summary csv is empty")

    current_problem = None
    for row in rows:
        if row["problem"] != current_problem:
            current_problem = row["problem"]
            print(f"\nproblem: {current_problem}")
            print("cfl\tfinest L1\tfinest L2\tfinest Linf\tlast order(L1)\tlast order(L2)\tlast order(Linf)")

        print(
            f"{float(row['cfl']):.3f}\t"
            f"{float(row['finest_l1']):.6e}\t"
            f"{float(row['finest_l2']):.6e}\t"
            f"{float(row['finest_linf']):.6e}\t"
            f"{row['last_order_l1'] or '-'}\t"
            f"{row['last_order_l2'] or '-'}\t"
            f"{row['last_order_linf'] or '-'}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
