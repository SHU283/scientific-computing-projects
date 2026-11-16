#!/usr/bin/env python3
"""
Read a convergence-summary CSV and print a compact summary.
If matplotlib is installed, the same data can also be plotted later.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "nx": float(row["nx"]),
                    "l1_rho": float(row["l1_rho"]),
                    "order_l1": float(row["order_l1"]),
                    "l2_rho": float(row["l2_rho"]),
                    "order_l2": float(row["order_l2"]),
                    "linf_rho": float(row["linf_rho"]),
                    "order_linf": float(row["order_linf"]),
                }
            )
    return rows


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python3 scripts/plot_convergence.py results/xxx_convergence.csv")
        return 1

    path = Path(sys.argv[1])
    if not path.exists():
        print(f"file not found: {path}")
        return 1

    rows = load_rows(path)
    print(f"Convergence summary for: {path}")
    for row in rows:
        nx = int(row["nx"])
        print(
            f"nx={nx:4d}  "
            f"L1={row['l1_rho']:.6e}  p1={row['order_l1']:.3f}  "
            f"L2={row['l2_rho']:.6e}  p2={row['order_l2']:.3f}  "
            f"Linf={row['linf_rho']:.6e}  pinf={row['order_linf']:.3f}"
        )

    if len(rows) >= 2:
        last = rows[-1]
        print("\nLatest refinement snapshot:")
        print(
            f"L1 order ~ {last['order_l1']:.3f}, "
            f"L2 order ~ {last['order_l2']:.3f}, "
            f"Linf order ~ {last['order_linf']:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

