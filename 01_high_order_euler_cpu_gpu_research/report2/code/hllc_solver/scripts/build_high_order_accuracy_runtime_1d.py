#!/usr/bin/env python3
"""
Merge 1D high-order convergence summaries with runtime summaries.

Current scope:
- problems: smooth_advection, gaussian_advection
- methods: muscl, weno2, weno3, ader3
- backends:
  - muscl, weno3, ader3: serial, openmp
  - weno2: serial, openmp, cuda
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


METHODS = ["muscl", "weno2", "weno3", "ader3"]
PROBLEMS = ["smooth_advection", "gaussian_advection"]


def load_convergence(path: Path) -> dict[int, dict[str, float]]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        data: dict[int, dict[str, float]] = {}
        for row in reader:
            nx = int(row["nx"])
            data[nx] = {
                "l1_rho": float(row["l1_rho"]),
                "l2_rho": float(row["l2_rho"]),
                "linf_rho": float(row["linf_rho"]),
            }
    return data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime", default="results/stageC/high_order_runtime_1d.csv")
    parser.add_argument("--output", default="results/stageC/high_order_accuracy_runtime_1d.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    runtime_path = root / args.runtime
    output_path = root / args.output

    if not runtime_path.exists():
        raise SystemExit(f"file not found: {runtime_path}")

    conv_by_key: dict[tuple[str, str], dict[int, dict[str, float]]] = {}
    for problem in PROBLEMS:
        for method in METHODS:
            conv_path = root / "results" / "stageA" / f"{problem}_1d_{method}_convergence.csv"
            if conv_path.exists():
                conv_by_key[(problem, method)] = load_convergence(conv_path)

    with runtime_path.open("r", encoding="utf-8") as f:
        runtime_rows = list(csv.DictReader(f))

    merged_rows: list[dict[str, str | int | float]] = []
    for row in runtime_rows:
        key = (row["problem"], row["method"])
        nx = int(row["nx"])
        conv = conv_by_key.get(key)
        if not conv or nx not in conv:
            continue

        merged_rows.append(
            {
                "problem": row["problem"],
                "backend": row["backend"],
                "method": row["method"],
                "nx": nx,
                "best_seconds": float(row["best_seconds"]),
                "mean_seconds": float(row["mean_seconds"]),
                "median_seconds": float(row["median_seconds"]),
                "l1_rho": conv[nx]["l1_rho"],
                "l2_rho": conv[nx]["l2_rho"],
                "linf_rho": conv[nx]["linf_rho"],
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "problem",
                "backend",
                "method",
                "nx",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
                "l1_rho",
                "l2_rho",
                "linf_rho",
            ],
        )
        writer.writeheader()
        writer.writerows(merged_rows)

    print(f"saved high-order accuracy-runtime csv: {output_path}")
    print(f"rows written: {len(merged_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
