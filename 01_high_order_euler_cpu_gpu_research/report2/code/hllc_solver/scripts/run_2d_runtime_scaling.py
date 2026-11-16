#!/usr/bin/env python3
"""Run the 2D Gaussian runtime-scaling benchmark.

Coverage:
- problem: gaussian_advection_2d
- methods: godunov, muscl
- backends: serial, openmp, cuda
- resolutions: 50, 100, 200

Output:
- results/stageC/runtime_scaling_2d.csv
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


METHODS = ["godunov", "muscl"]
BACKENDS = ["serial", "openmp", "cuda"]
RESOLUTIONS = [50, 100, 200]


def pick_executable(backend: str) -> str:
    return "./build-cuda/final_project" if backend == "cuda" else "./build/final_project"


def timed_run(cmd: list[str], cwd: Path) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    end = time.perf_counter()
    return end - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", default="results/stageC/runtime_scaling_2d.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | float | int]] = []

    for n in RESOLUTIONS:
        for method in METHODS:
            for backend in BACKENDS:
                exe = pick_executable(backend)
                output_csv = f"results/stageC/runtime_scaling_2d_{backend}_{method}_n{n}.csv"
                cmd = [
                    exe,
                    f"--backend={backend}",
                    "--problem=gaussian_advection_2d",
                    f"--method={method}",
                    f"--nx={n}",
                    f"--ny={n}",
                    "--final_time=0.1",
                    f"--output={output_csv}",
                ]

                print(f"\n[scaling] n={n} backend={backend} method={method}")
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
                        "problem": "gaussian_advection_2d",
                        "backend": backend,
                        "method": method,
                        "nx": n,
                        "ny": n,
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
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nsaved runtime scaling csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

