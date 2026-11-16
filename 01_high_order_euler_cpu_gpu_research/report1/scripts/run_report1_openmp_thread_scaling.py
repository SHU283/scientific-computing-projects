#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ["muscl", "weno2"]
THREAD_COUNTS = [1, 2, 4, 8, 16]


def find_executable() -> Path:
    candidates = [
        ROOT / "build-cuda" / "final_project",
        ROOT / "build" / "final_project",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit("No final_project executable found in build-cuda/ or build/.")


def timed_run(cmd: list[str], threads: int) -> float:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("OMP_PROC_BIND", "false")

    start = time.perf_counter()
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env)
    end = time.perf_counter()
    return end - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=100)
    parser.add_argument("--ny", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", default="results/stageC/openmp_thread_scaling_report1.csv")
    args = parser.parse_args()

    executable = find_executable()
    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []
    for method in METHODS:
        for threads in THREAD_COUNTS:
            output_csv = (
                f"results/stageC/bench_openmp_thread_scaling_"
                f"gaussian_advection_2d_{method}_n{args.nx}_t{threads}.csv"
            )
            cmd = [
                str(executable),
                "--backend=openmp",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={args.nx}",
                f"--ny={args.ny}",
                "--boundary=periodic",
                "--final_time=0.1",
                f"--output={output_csv}",
            ]

            print(f"\n[report1-openmp-thread-scaling] method={method} n={args.nx} threads={threads}")
            for warmup_idx in range(args.warmup):
                elapsed = timed_run(cmd, threads)
                print(f"  warmup {warmup_idx + 1}: {elapsed:.6f} s")

            samples: list[float] = []
            for repeat_idx in range(args.repeats):
                elapsed = timed_run(cmd, threads)
                samples.append(elapsed)
                print(f"  repeat {repeat_idx + 1}: {elapsed:.6f} s")

            rows.append(
                {
                    "problem": "gaussian_advection_2d",
                    "backend": "openmp",
                    "method": method,
                    "nx": args.nx,
                    "ny": args.ny,
                    "threads": threads,
                    "warmup": args.warmup,
                    "repeats": args.repeats,
                    "best_seconds": min(samples),
                    "mean_seconds": statistics.fmean(samples),
                    "median_seconds": statistics.median(samples),
                }
            )

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "problem",
                "backend",
                "method",
                "nx",
                "ny",
                "threads",
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nsaved Report 1 OpenMP thread-scaling csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

