#!/usr/bin/env python3
"""
Run 2D runtime-scaling benchmarks across serial, OpenMP, and CUDA backends.

Each case is repeated several times and summarized by best/mean/median runtime so
the resulting CSV can support backend speedup comparisons.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ["godunov", "muscl", "weno2"]
BACKENDS = ["serial", "openmp", "cuda"]
RESOLUTIONS = [50, 100, 200]


def find_executable() -> Path:
    candidates = [
        ROOT / "build-cuda" / "final_project",
        ROOT / "build" / "final_project",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit("No final_project executable found in build-cuda/ or build/.")


def timed_run(cmd: list[str]) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    end = time.perf_counter()
    return end - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", default="results/stageC/runtime_scaling_2d_report1_all_backends.csv")
    args = parser.parse_args()

    executable = find_executable()
    if "build-cuda" not in executable.as_posix():
        raise SystemExit("Please use a CUDA build first: cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON")

    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []
    for method in METHODS:
        for backend in BACKENDS:
            for n in RESOLUTIONS:
                output_csv = f"results/stageC/bench_gaussian_advection_2d_{backend}_{method}_n{n}.csv"
                cmd = [
                    str(executable),
                    f"--backend={backend}",
                    "--problem=gaussian_advection_2d",
                    f"--method={method}",
                    f"--nx={n}",
                    f"--ny={n}",
                    "--boundary=periodic",
                    "--final_time=0.1",
                    f"--output={output_csv}",
                ]

                print(f"\n[report1-2d-runtime] method={method} backend={backend} n={n}")
                for warmup_idx in range(args.warmup):
                    elapsed = timed_run(cmd)
                    print(f"  warmup {warmup_idx + 1}: {elapsed:.6f} s")

                samples: list[float] = []
                for repeat_idx in range(args.repeats):
                    elapsed = timed_run(cmd)
                    samples.append(elapsed)
                    print(f"  repeat {repeat_idx + 1}: {elapsed:.6f} s")

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

    print(f"\nsaved Report 1 all-backend runtime csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
