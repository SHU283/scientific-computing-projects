#!/usr/bin/env python3
"""Benchmark Report 2 OpenMP scaling on the smooth 2D Gaussian problem."""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("weno2", "weno3", "ader3")
THREAD_COUNTS = (1, 2, 4, 8, 16)
METHOD_CFL = {
    "weno2": 0.4,
    "weno3": 0.4,
    "ader3": 0.25,
}
FIELDNAMES = (
    "problem",
    "backend",
    "method",
    "nx",
    "ny",
    "cfl",
    "final_time",
    "threads",
    "warmup",
    "repeats",
    "best_seconds",
    "mean_seconds",
    "median_seconds",
    "stdev_seconds",
    "relative_stdev",
    "samples_seconds",
)


def find_executable(build_dir: str | None) -> Path:
    candidates = []
    if build_dir:
        candidates.append(ROOT / build_dir / "final_project")
    candidates.extend(
        [
            ROOT / "build-report2-cuda" / "final_project",
            ROOT / "build-cuda" / "final_project",
            ROOT / "build" / "final_project",
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit("No final_project executable found in a known build directory.")


def completed_keys(output_path: Path) -> set[tuple[str, int]]:
    if not output_path.exists():
        return set()
    with output_path.open(newline="", encoding="utf-8") as stream:
        return {
            (row["method"], int(row["threads"]))
            for row in csv.DictReader(stream)
        }


def timed_run(command: list[str], threads: int) -> float:
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(threads)
    environment["OMP_DYNAMIC"] = "FALSE"
    environment.setdefault("OMP_PROC_BIND", "close")
    environment.setdefault("OMP_PLACES", "cores")

    start = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        print(completed.stdout)
        raise subprocess.CalledProcessError(
            completed.returncode,
            command,
            output=completed.stdout,
        )
    return elapsed


def append_row(output_path: Path, row: dict[str, str | int | float]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not output_path.exists() or output_path.stat().st_size == 0
    with output_path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default=None)
    parser.add_argument("--nx", type=int, default=100)
    parser.add_argument("--ny", type=int, default=100)
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, nargs="+", default=THREAD_COUNTS)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument(
        "--output",
        default=(
            "results/report2/performance/openmp_thread_scaling/"
            "openmp_thread_scaling_report2_lovelace_a30.csv"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    executable = find_executable(args.build_dir)
    output_path = ROOT / args.output
    done = completed_keys(output_path) if args.resume else set()

    for method in args.methods:
        cfl = METHOD_CFL[method]
        for threads in args.threads:
            key = (method, threads)
            if key in done:
                print(
                    f"\n[report2-openmp-thread-scaling] skipping completed "
                    f"method={method} threads={threads}"
                )
                continue

            command = [
                str(executable),
                "--backend=openmp",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={args.nx}",
                f"--ny={args.ny}",
                "--boundary=periodic",
                f"--final_time={args.final_time}",
                f"--cfl={cfl}",
                "--output=none",
            ]

            print(
                f"\n[report2-openmp-thread-scaling] method={method} "
                f"n={args.nx} threads={threads} cfl={cfl}"
            )
            for warmup_index in range(args.warmup):
                elapsed = timed_run(command, threads)
                print(f"  warmup {warmup_index + 1}: {elapsed:.6f} s")

            samples: list[float] = []
            for repeat_index in range(args.repeats):
                elapsed = timed_run(command, threads)
                samples.append(elapsed)
                print(f"  repeat {repeat_index + 1}: {elapsed:.6f} s")

            mean_seconds = statistics.fmean(samples)
            stdev_seconds = statistics.stdev(samples) if len(samples) > 1 else 0.0
            append_row(
                output_path,
                {
                    "problem": "gaussian_advection_2d",
                    "backend": "openmp",
                    "method": method,
                    "nx": args.nx,
                    "ny": args.ny,
                    "cfl": cfl,
                    "final_time": args.final_time,
                    "threads": threads,
                    "warmup": args.warmup,
                    "repeats": args.repeats,
                    "best_seconds": min(samples),
                    "mean_seconds": mean_seconds,
                    "median_seconds": statistics.median(samples),
                    "stdev_seconds": stdev_seconds,
                    "relative_stdev": (
                        stdev_seconds / mean_seconds if mean_seconds else 0.0
                    ),
                    "samples_seconds": ";".join(
                        f"{sample:.9g}" for sample in samples
                    ),
                },
            )
            print(f"  saved row to: {output_path}")

    print(f"\nsaved Report 2 OpenMP thread-scaling csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
