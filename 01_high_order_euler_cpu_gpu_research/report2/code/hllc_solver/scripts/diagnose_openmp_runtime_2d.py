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
METHODS = ["godunov", "muscl", "weno2"]


def find_executable() -> Path:
    for candidate in [ROOT / "build-cuda" / "final_project", ROOT / "build" / "final_project"]:
        if candidate.exists():
            return candidate
    raise SystemExit("No final_project executable found.")


def timed_run(cmd: list[str], env: dict[str, str]) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, env=env)
    return time.perf_counter() - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", nargs="+", default=["default", "1", "2", "4", "8", "16"])
    parser.add_argument("--resolutions", type=int, nargs="+", default=[25, 50, 100])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--output", default="results/stageC/rerun_2026-05-24/openmp_thread_diagnostic_25_50_100.csv")
    args = parser.parse_args()

    executable = find_executable()
    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []
    for thread_label in args.threads:
        for method in METHODS:
            for n in args.resolutions:
                env = os.environ.copy()
                env["OMP_DYNAMIC"] = "FALSE"
                if thread_label == "default":
                    env.pop("OMP_NUM_THREADS", None)
                else:
                    env["OMP_NUM_THREADS"] = thread_label

                output_csv = (
                    f"results/stageC/rerun_2026-05-24/"
                    f"diag_openmp_{thread_label}_{method}_n{n}.csv"
                )
                cmd = [
                    str(executable),
                    "--backend=openmp",
                    "--problem=gaussian_advection_2d",
                    f"--method={method}",
                    f"--nx={n}",
                    f"--ny={n}",
                    "--boundary=periodic",
                    "--final_time=0.1",
                    f"--output={output_csv}",
                ]

                print(f"\n[openmp-diagnostic] threads={thread_label} method={method} n={n}")
                for warmup_idx in range(args.warmup):
                    elapsed = timed_run(cmd, env)
                    print(f"  warmup {warmup_idx + 1}: {elapsed:.6f} s")

                samples: list[float] = []
                for repeat_idx in range(args.repeats):
                    elapsed = timed_run(cmd, env)
                    samples.append(elapsed)
                    print(f"  repeat {repeat_idx + 1}: {elapsed:.6f} s")

                rows.append(
                    {
                        "backend": "openmp",
                        "threads": thread_label,
                        "method": method,
                        "nx": n,
                        "ny": n,
                        "warmup": args.warmup,
                        "repeats": args.repeats,
                        "best_seconds": min(samples),
                        "mean_seconds": statistics.fmean(samples),
                        "median_seconds": statistics.median(samples),
                        "samples": ";".join(f"{sample:.9f}" for sample in samples),
                    }
                )

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "backend",
                "threads",
                "method",
                "nx",
                "ny",
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
                "samples",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nsaved diagnostic csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
