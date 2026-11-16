#!/usr/bin/env python3
"""
Run Report 1 WENO-2 2D discontinuous artifacts:
- Stage A: representative single runs and field plots
- Stage B: serial/openmp consistency
- Stage C: runtime benchmarks for discontinuous 2D cases
- Stage D: three-method comparison plots with godunov/muscl/weno2
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


CASES = [
    {
        "label": "quadrant_riemann_2d",
        "problem": "quadrant_riemann_2d",
        "nx": 100,
        "ny": 100,
        "final_time": 0.2,
        "field": "rho",
    },
    {
        "label": "shock_bubble_2d",
        "problem": "shock_bubble_2d",
        "nx": 100,
        "ny": 100,
        "final_time": 0.05,
        "field": "rho",
    },
]

BACKENDS = ["serial", "openmp", "cuda"]


def pick_executable(backend: str) -> str:
    if backend == "cuda":
        return "./build-cuda/final_project"
    return "./build/final_project"


def run(cmd: list[str], cwd: Path) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)
    print()


def timed_run(cmd: list[str], cwd: Path) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return time.perf_counter() - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--runtime-output", default="results/stageC/runtime_discontinuous_2d_weno2.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    rows: list[dict[str, str | int | float]] = []

    for case in CASES:
        label = case["label"]
        problem = case["problem"]
        nx = case["nx"]
        ny = case["ny"]
        final_time = case["final_time"]
        field = case["field"]

        # Stage A representative run.
        stage_a_csv = f"results/stageA/{label}_weno2.csv"
        run(
            [
                pick_executable("serial"),
                "--backend=serial",
                f"--problem={problem}",
                "--method=weno2",
                f"--nx={nx}",
                f"--ny={ny}",
                "--boundary=outflow",
                f"--final_time={final_time}",
                f"--output={stage_a_csv}",
            ],
            root,
        )
        run(["python3", "scripts/plot_field_2d.py", stage_a_csv, f"--field={field}"], root)

        # Stage B consistency.
        serial_csv = f"results/stageB/{label}_serial_weno2.csv"
        openmp_csv = f"results/stageB/{label}_openmp_weno2.csv"
        for backend, output_csv in [("serial", serial_csv), ("openmp", openmp_csv)]:
            run(
                [
                    pick_executable(backend),
                    f"--backend={backend}",
                    f"--problem={problem}",
                    "--method=weno2",
                    f"--nx={nx}",
                    f"--ny={ny}",
                    "--boundary=outflow",
                    f"--final_time={final_time}",
                    f"--output={output_csv}",
                ],
                root,
            )
        run(["python3", "scripts/compare_csv_results.py", serial_csv, openmp_csv], root)

        # Stage D comparison vs baseline methods.
        run(
            [
                "python3",
                "scripts/compare_fields_2d_multi.py",
                f"results/stageA/{label}_godunov.csv",
                f"results/stageA/{label}_muscl.csv",
                stage_a_csv,
                "--field",
                field,
            ],
            root,
        )

        # Stage C runtime.
        for backend in BACKENDS:
            output_csv = f"results/stageC/bench_{label}_{backend}_weno2.csv"
            cmd = [
                pick_executable(backend),
                f"--backend={backend}",
                f"--problem={problem}",
                "--method=weno2",
                f"--nx={nx}",
                f"--ny={ny}",
                f"--final_time={final_time}",
                f"--output={output_csv}",
            ]
            print(f"[weno2-2d-discontinuous] case={label} backend={backend}")
            for warmup in range(args.warmup):
                elapsed = timed_run(cmd, root)
                print(f"  warmup {warmup + 1}: {elapsed:.6f} s")
            samples: list[float] = []
            for repeat in range(args.repeats):
                elapsed = timed_run(cmd, root)
                samples.append(elapsed)
                print(f"  repeat {repeat + 1}: {elapsed:.6f} s")
            rows.append(
                {
                    "case": label,
                    "problem": problem,
                    "backend": backend,
                    "method": "weno2",
                    "nx": nx,
                    "ny": ny,
                    "final_time": final_time,
                    "warmup": args.warmup,
                    "repeats": args.repeats,
                    "best_seconds": min(samples),
                    "mean_seconds": statistics.fmean(samples),
                    "median_seconds": statistics.median(samples),
                }
            )
            print()

    output_path = root / args.runtime_output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "problem",
                "backend",
                "method",
                "nx",
                "ny",
                "final_time",
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"saved WENO-2 discontinuous runtime csv: {output_path}")
    print("WENO-2 2D discontinuous suite finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

