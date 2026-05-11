#!/usr/bin/env python3
"""
Build Report 1 accuracy-runtime table for 2D smooth methods:
godunov, muscl, weno2.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_convergence(path: Path) -> dict[int, dict[str, float]]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return {
            int(row["n"]): {
                "l1_rho": float(row["l1_rho"]),
                "l2_rho": float(row["l2_rho"]),
                "linf_rho": float(row["linf_rho"]),
            }
            for row in reader
        }


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime", default="results/stageC/runtime_scaling_2d_report1_all_backends.csv")
    parser.add_argument("--godunov", default="results/stageA/gaussian_advection_2d_godunov_convergence.csv")
    parser.add_argument("--muscl", default="results/stageA/gaussian_advection_2d_muscl_convergence.csv")
    parser.add_argument("--weno2", default="results/stageA/gaussian_advection_2d_weno2_convergence.csv")
    parser.add_argument("--output", default="results/stageC/accuracy_runtime_2d_report1.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    conv_by_method = {
        "godunov": load_convergence(root / args.godunov),
        "muscl": load_convergence(root / args.muscl),
        "weno2": load_convergence(root / args.weno2),
    }
    runtime_rows = load_rows(root / args.runtime)

    merged_rows: list[dict[str, str | int | float]] = []
    for row in runtime_rows:
        method = row["method"]
        n = int(row["nx"])
        conv = conv_by_method.get(method, {})
        if n not in conv:
            continue
        merged_rows.append(
            {
                "problem": row["problem"],
                "backend": row["backend"],
                "method": method,
                "nx": n,
                "ny": int(row["ny"]),
                "best_seconds": float(row["best_seconds"]),
                "mean_seconds": float(row["mean_seconds"]),
                "median_seconds": float(row["median_seconds"]),
                "l1_rho": conv[n]["l1_rho"],
                "l2_rho": conv[n]["l2_rho"],
                "linf_rho": conv[n]["linf_rho"],
            }
        )

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

    print(f"saved Report 1 accuracy-runtime csv: {output_path}")
    print(f"rows written: {len(merged_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
