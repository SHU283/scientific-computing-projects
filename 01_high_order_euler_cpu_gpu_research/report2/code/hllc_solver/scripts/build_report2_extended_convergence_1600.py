#!/usr/bin/env python3
"""Merge the canonical 1D sequence with the verified n=1600 extension."""

from __future__ import annotations

import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    canonical = (
        root
        / "results/report2_hllc/cfl_aligned_validation/validation_1d/"
        "report2_high_order_convergence_1d.csv"
    )
    extension = (
        root
        / "results/report2_hllc/extended_convergence_1600/"
        "validation_1d_openmp/report2_high_order_convergence_1d.csv"
    )
    output = (
        root
        / "results/report2_hllc/extended_convergence_1600/"
        "validation_1d_openmp/"
        "report2_high_order_convergence_1d_50_1600.csv"
    )

    rows_by_key: dict[tuple[str, str, int], dict[str, str]] = {}
    for row in read_rows(canonical) + read_rows(extension):
        key = (row["problem"], row["method"], int(row["nx"]))
        rows_by_key[key] = row

    rows = [
        rows_by_key[key]
        for key in sorted(rows_by_key, key=lambda item: (item[0], item[1], item[2]))
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
