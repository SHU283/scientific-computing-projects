#!/usr/bin/env python3
"""
Run Stage C runtime benchmarks for 2D discontinuous problems.

Current coverage:
- quadrant_riemann_2d
- shock_bubble_2d

Methods:
- godunov
- muscl
- weno2

Backends:
- serial
- openmp
- cuda
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
        "label": "quadrant_riemann_2d",
        "problem": "quadrant_riemann_2d",
        "nx": 100,
        "ny": 100,
        "final_time": 0.2,
    },
    {
        "label": "shock_bubble_2d",
        "problem": "shock_bubble_2d",
        "nx": 100,
        "ny": 100,
        "final_time": 0.05,
    },
]

METHODS = ["godunov", "muscl", "weno2"]
BACKENDS = ["serial", "openmp", "cuda"]


def pick_executable(root: Path, backend: str) -> str:
    if backend == "cuda":
        return "./build-cuda/final_project"
    return "./build/final_project"


def timed_run(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env)
    end = time.perf_counter()
    return end - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--nx", type=int, default=None, help="Override the x-resolution for all discontinuous cases.")
    parser.add_argument("--ny", type=int, default=None, help="Override the y-resolution for all discontinuous cases.")
    parser.add_argument("--openmp-threads", type=int, default=16)
    parser.add_argument("--write-fields", action="store_true", help="Write field CSV files during timed runs. By default timing uses --output=none.")
    parser.add_argument("--output", default="results/stageC/runtime_discontinuous_2d_report1_all_backends_noio.csv")
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
        nx = args.nx if args.nx is not None else case["nx"]
        ny = args.ny if args.ny is not None else case["ny"]
        for method in METHODS:
            for backend in BACKENDS:
                exe = pick_executable(root, backend)
                output_csv = f"results/stageC/bench_{case['label']}_{backend}_{method}.csv" if args.write_fields else "none"
                cmd = [
                    exe,
                    f"--backend={backend}",
                    f"--problem={case['problem']}",
                    f"--method={method}",
                    f"--nx={nx}",
                    f"--ny={ny}",
                    f"--final_time={case['final_time']}",
                    f"--output={output_csv}",
                ]

                print(f"\n[benchmark] case={case['label']} backend={backend} method={method}")
                run_env = None
                if backend == "openmp":
                    import os

                    run_env = os.environ.copy()
                    run_env["OMP_NUM_THREADS"] = str(args.openmp_threads)
                    run_env["OMP_DYNAMIC"] = "FALSE"

                for warmup in range(args.warmup):
                    elapsed = timed_run(cmd, root, run_env)
                    print(f"  warmup {warmup + 1}: {elapsed:.6f} s")

                samples: list[float] = []
                for repeat in range(args.repeats):
                    elapsed = timed_run(cmd, root, run_env)
                    samples.append(elapsed)
                    print(f"  repeat {repeat + 1}: {elapsed:.6f} s")

                rows.append(
                    {
                        "case": case["label"],
                        "problem": case["problem"],
                        "backend": backend,
                        "method": method,
                        "nx": nx,
                        "ny": ny,
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

    print(f"\nsaved discontinuous runtime csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
