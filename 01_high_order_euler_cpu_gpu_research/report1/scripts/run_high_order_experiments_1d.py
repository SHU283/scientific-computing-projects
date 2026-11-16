#!/usr/bin/env python3
"""
Run a compact 1D experimental high-order comparison on gaussian_advection.

Current methods:
- muscl
- weno2
- weno3

Outputs:
- solver-generated convergence CSVs in results/stageA
- a compact summary CSV in results/stageA
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


METHODS = ["muscl", "weno2", "weno3"]


def parse_last_observed_order(csv_path: Path) -> tuple[float | None, float | None, float | None]:
    with csv_path.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if len(rows) < 2:
        return None, None, None

    last = rows[-1]

    def maybe_float(key: str) -> float | None:
        value = last.get(key, "")
        if value in ("", "-", None):
            return None
        return float(value)

    return (
        maybe_float("order_l1"),
        maybe_float("order_l2"),
        maybe_float("order_linf"),
    )


def parse_finest_error(csv_path: Path) -> tuple[float, float, float]:
    with csv_path.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise RuntimeError(f"empty convergence csv: {csv_path}")

    last = rows[-1]
    return float(last["l1_rho"]), float(last["l2_rho"]), float(last["linf_rho"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--problem", default="gaussian_advection", choices=["gaussian_advection", "smooth_advection"])
    parser.add_argument("--output", default="results/stageA/high_order_experiments_1d.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | float]] = []

    for method in METHODS:
        cmd = [
            "./build/final_project",
            "--study=convergence",
            f"--problem={args.problem}",
            f"--method={method}",
        ]

        print(f"\n[high-order] problem={args.problem} method={method}")
        subprocess.run(cmd, cwd=root, check=True)

        convergence_csv = root / "results" / "stageA" / f"{args.problem}_1d_{method}_convergence.csv"
        l1, l2, linf = parse_finest_error(convergence_csv)
        order_l1, order_l2, order_linf = parse_last_observed_order(convergence_csv)

        rows.append(
            {
                "problem": args.problem,
                "method": method,
                "finest_l1": l1,
                "finest_l2": l2,
                "finest_linf": linf,
                "last_order_l1": "" if order_l1 is None else order_l1,
                "last_order_l2": "" if order_l2 is None else order_l2,
                "last_order_linf": "" if order_linf is None else order_linf,
            }
        )

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "problem",
                "method",
                "finest_l1",
                "finest_l2",
                "finest_linf",
                "last_order_l1",
                "last_order_l2",
                "last_order_linf",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nsaved high-order summary csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

