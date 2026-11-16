#!/usr/bin/env python3
"""Merge benchmark CSV files, replacing duplicate configuration rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


KEY_FIELDS = ("problem", "method", "backend", "nx", "ny", "cfl", "final_time")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def row_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in KEY_FIELDS)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_csv")
    parser.add_argument("replacement_csv")
    parser.add_argument("output_csv")
    args = parser.parse_args()

    base_path = Path(args.base_csv)
    replacement_path = Path(args.replacement_csv)
    output_path = Path(args.output_csv)
    base_rows = read_rows(base_path)
    replacement_rows = read_rows(replacement_path)

    merged = {row_key(row): row for row in base_rows}
    for row in replacement_rows:
        merged[row_key(row)] = row

    rows = sorted(
        merged.values(),
        key=lambda row: (
            row["problem"],
            row["method"],
            row["backend"],
            int(row["nx"]),
        ),
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"merged {len(base_rows)} base rows and "
        f"{len(replacement_rows)} replacement rows into {len(rows)} rows"
    )
    print(f"saved merged CSV: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
