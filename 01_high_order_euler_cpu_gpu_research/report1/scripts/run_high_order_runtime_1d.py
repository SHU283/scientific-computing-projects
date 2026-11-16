#!/usr/bin/env python3
"""
Benchmark current high-order 1D methods on supported CPU backends.

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
import statistics
import subprocess
import time
from pathlib import Path


CASES = [
    {
        "label": "1d_smooth",
        "problem": "smooth_advection",
        "nx": 400,
        "final_time": 0.3,
    },
    {
        "label": "1d_gaussian",
        "problem": "gaussian_advection",
        "nx": 400,
        "final_time": 0.3,
    },
]

METHODS = ["muscl", "weno2", "weno3", "ader3"]


def backends_for_method(method: str) -> list[str]:
    if method == "weno2":
        return ["serial", "openmp", "cuda"]
    return ["serial", "openmp"]


def executable_for_backend(backend: str) -> str:
    if backend == "cuda":
        return "./build-cuda/final_project"
    return "./build/final_project"


def timed_run(cmd: list[str], cwd: Path) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    end = time.perf_counter()
    return end - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", default="results/stageC/high_order_runtime_1d.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []

    for case in CASES:
        for method in METHODS:
            for backend in backends_for_method(method):
                output_csv = f"results/stageC/bench_{case['label']}_{backend}_{method}.csv"
                cmd = [
                    executable_for_backend(backend),
                    f"--backend={backend}",
                    f"--problem={case['problem']}",
                    f"--method={method}",
                    f"--nx={case['nx']}",
                    "--boundary=periodic",
                    f"--final_time={case['final_time']}",
                    f"--output={output_csv}",
                ]

                print(f"\n[high-order-benchmark] case={case['label']} backend={backend} method={method}")
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
                        "case": case["label"],
                        "problem": case["problem"],
                        "backend": backend,
                        "method": method,
                        "nx": case["nx"],
                        "final_time": case["final_time"],
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
                "case",
                "problem",
                "backend",
                "method",
                "nx",
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

    print(f"\nsaved high-order runtime csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

