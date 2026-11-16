#!/usr/bin/env python3
"""
Build Report 1 discontinuous 2D runtime tables from baseline and WENO-2 sources.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BASE_CSV = ROOT / "results" / "stageC" / "runtime_discontinuous_2d.csv"
WENO2_ALIGNMENT_CSV = ROOT / "results" / "stageC" / "runtime_weno2_backend_alignment.csv"
WENO2_OUTPUT_CSV = ROOT / "results" / "stageC" / "runtime_discontinuous_2d_weno2.csv"
OUTPUT_CSV = ROOT / "results" / "stageC" / "runtime_discontinuous_2d_report1_all_backends.csv"


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"missing csv: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def build_weno2_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    valid_cases = {"quadrant_riemann_2d", "shock_bubble_2d"}
    for row in rows:
        problem = row["problem"]
        if problem not in valid_cases:
            continue
        selected.append(
            {
                "case": problem,
                "problem": problem,
                "backend": row["backend"],
                "method": row["method"],
                "nx": row["nx"],
                "ny": row["ny"],
                "final_time": row["final_time"],
                "warmup": row["warmup"],
                "repeats": row["repeats"],
                "best_seconds": row["best_seconds"],
                "mean_seconds": row["mean_seconds"],
                "median_seconds": row["median_seconds"],
            }
        )
    return selected


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    base_rows = load_rows(BASE_CSV)
    weno2_alignment_rows = load_rows(WENO2_ALIGNMENT_CSV)
    weno2_rows = build_weno2_rows(weno2_alignment_rows)
    rows = base_rows + weno2_rows
    if not rows:
        raise SystemExit("no rows found to merge")

    fieldnames = list(rows[0].keys())
    write_rows(WENO2_OUTPUT_CSV, weno2_rows, fieldnames)
    write_rows(OUTPUT_CSV, rows, fieldnames)

    print(f"saved rebuilt WENO-2 discontinuous runtime csv: {WENO2_OUTPUT_CSV}")
    print(f"rows written: {len(weno2_rows)}")
    print(f"saved merged discontinuous runtime csv: {OUTPUT_CSV}")
    print(f"rows written: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

