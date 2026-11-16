#!/usr/bin/env python3
"""Summarize timing variability and backend speedups from a benchmark CSV."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


KEY_FIELDS = ("problem", "method", "nx")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path")
    parser.add_argument(
        "--metric",
        choices=("best_seconds", "median_seconds"),
        default="median_seconds",
    )
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    rows = read_rows(Path(args.csv_path))
    print("Highest relative standard deviations")
    print("problem,method,backend,n,repeats,relative_stdev")
    for row in sorted(
        rows,
        key=lambda item: float(item["relative_stdev"]),
        reverse=True,
    )[: args.top]:
        print(
            f"{row['problem']},{row['method']},{row['backend']},"
            f"{row['nx']},{row['repeats']},"
            f"{float(row['relative_stdev']):.3f}"
        )

    grouped: dict[tuple[str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for row in rows:
        problem, method, resolution = (row[field] for field in KEY_FIELDS)
        grouped[(problem, method, resolution)][row["backend"]] = row

    print(f"\nBackend summary using {args.metric}")
    print(
        "problem,method,n,serial_seconds,openmp_seconds,cuda_seconds,"
        "openmp_speedup,cuda_speedup"
    )
    for key in sorted(grouped, key=lambda item: (item[0], item[1], int(item[2]))):
        backends = grouped[key]
        if "serial" not in backends:
            continue
        serial = float(backends["serial"][args.metric])
        openmp = (
            float(backends["openmp"][args.metric])
            if "openmp" in backends
            else float("nan")
        )
        cuda = (
            float(backends["cuda"][args.metric])
            if "cuda" in backends
            else float("nan")
        )
        print(
            f"{key[0]},{key[1]},{key[2]},"
            f"{serial:.6g},{openmp:.6g},{cuda:.6g},"
            f"{serial / openmp:.3f},{serial / cuda:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
