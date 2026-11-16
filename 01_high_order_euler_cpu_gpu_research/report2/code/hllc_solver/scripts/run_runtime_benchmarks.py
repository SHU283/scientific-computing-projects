#!/usr/bin/env python3
"""Run the baseline runtime campaign and export its measurements to CSV.

Coverage:
- problems: 1D smooth_advection and 2D gaussian_advection_2d
- methods: godunov and muscl
- backends: serial, openmp, and cuda

Each configuration receives one warm-up and five measured repetitions. Timing
uses external wall-clock measurements so that every backend is compared through
the same boundary. The CSV records minimum, mean, and median runtimes.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


BENCH_CASES = [
    {
        "label": "1d_smooth",
        "problem": "smooth_advection",
        "nx": 400,
        "ny": 1,
        "final_time": 0.3,
    },
    {
        "label": "2d_gaussian",
        "problem": "gaussian_advection_2d",
        "nx": 100,
        "ny": 100,
        "final_time": 0.1,
    },
]

METHODS = ["godunov", "muscl"]
BACKENDS = ["serial", "openmp", "cuda"]


def pick_executable(root: Path, backend: str) -> str:
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
    parser.add_argument("--output", default="results/stageC/runtime_benchmarks.csv")
    args = parser.parse_args()

    if args.repeats <= 0:
        raise SystemExit("repeats must be positive")
    if args.warmup < 0:
        raise SystemExit("warmup must be non-negative")

    root = Path(__file__).resolve().parent.parent
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | float | int]] = []

    for case in BENCH_CASES:
        for method in METHODS:
            for backend in BACKENDS:
                exe = pick_executable(root, backend)
                output_csv = (
                    f"results/stageC/bench_{case['label']}_{backend}_{method}.csv"
                )
                cmd = [
                    exe,
                    f"--backend={backend}",
                    f"--problem={case['problem']}",
                    f"--method={method}",
                    f"--nx={case['nx']}",
                    f"--ny={case['ny']}",
                    f"--final_time={case['final_time']}",
                    f"--output={output_csv}",
                ]

                print(f"\n[benchmark] case={case['label']} backend={backend} method={method}")
                for warmup in range(args.warmup):
                    elapsed = timed_run(cmd, root)
                    print(f"  warmup {warmup + 1}: {elapsed:.6f} s")

                samples: list[float] = []
                for repeat in range(args.repeats):
                    elapsed = timed_run(cmd, root)
                    samples.append(elapsed)
                    print(f"  repeat {repeat + 1}: {elapsed:.6f} s")

                best = min(samples)
                mean = statistics.fmean(samples)
                median = statistics.median(samples)
                rows.append(
                    {
                        "case": case["label"],
                        "problem": case["problem"],
                        "backend": backend,
                        "method": method,
                        "nx": case["nx"],
                        "ny": case["ny"],
                        "final_time": case["final_time"],
                        "warmup": args.warmup,
                        "repeats": args.repeats,
                        "best_seconds": best,
                        "mean_seconds": mean,
                        "median_seconds": median,
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

    print(f"\nsaved benchmark csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

