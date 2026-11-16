#!/usr/bin/env python3
"""Benchmark Report 2 methods on the two discontinuous 2D problems."""

from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("weno2", "weno3", "ader3")
BACKENDS = ("serial", "openmp", "cuda")
CASES = {
    "quadrant_riemann_2d": {
        "boundary": "outflow",
        "final_time": 0.2,
    },
    "shock_bubble_2d": {
        "boundary": "outflow",
        "final_time": 0.2,
    },
}
METHOD_CFL = {
    "weno2": 0.2,
    "weno3": 0.2,
    "ader3": 0.14,
}
STEP_PATTERN = re.compile(r"time steps\s+:\s+(?P<steps>\d+)")
FIELDNAMES = (
    "problem",
    "method",
    "backend",
    "nx",
    "ny",
    "cfl",
    "final_time",
    "openmp_threads",
    "steps",
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
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit(
        "No CUDA-enabled final_project executable found. "
        "Build with -DFP_ENABLE_CUDA=ON first."
    )


def timed_run(
    command: list[str],
    environment: dict[str, str],
) -> tuple[float, str]:
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
    return elapsed, completed.stdout


def parse_steps(output: str) -> int:
    match = STEP_PATTERN.search(output)
    if not match:
        raise RuntimeError(f"could not parse time-step count:\n{output}")
    return int(match.group("steps"))


def completed_keys(output_path: Path) -> set[tuple[str, str, str, int]]:
    if not output_path.exists():
        return set()
    with output_path.open(newline="", encoding="utf-8") as stream:
        return {
            (
                row["problem"],
                row["method"],
                row["backend"],
                int(row["nx"]),
            )
            for row in csv.DictReader(stream)
        }


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
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--resolutions",
        type=int,
        nargs="+",
        default=(25, 50, 100, 200, 600),
    )
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--backends", nargs="+", choices=BACKENDS, default=BACKENDS)
    parser.add_argument("--cases", nargs="+", choices=tuple(CASES), default=tuple(CASES))
    parser.add_argument("--openmp-threads", type=int, default=16)
    parser.add_argument("--quadrant-final-time", type=float, default=0.2)
    parser.add_argument("--shock-final-time", type=float, default=0.2)
    parser.add_argument(
        "--output",
        default=(
            "results/report2/performance/"
            "runtime_discontinuous_2d_report2_lovelace_a30.csv"
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="skip rows that are already present in the output CSV",
    )
    args = parser.parse_args()

    if args.warmup < 0:
        raise SystemExit("warmup must be non-negative")
    if args.repeats <= 0:
        raise SystemExit("repeats must be positive")
    if any(resolution <= 0 for resolution in args.resolutions):
        raise SystemExit("resolutions must be positive")

    executable = find_executable(args.build_dir)
    output_path = ROOT / args.output
    done = completed_keys(output_path) if args.resume else set()
    final_times = {
        "quadrant_riemann_2d": args.quadrant_final_time,
        "shock_bubble_2d": args.shock_final_time,
    }

    for case_name in args.cases:
        case = CASES[case_name]
        for method in args.methods:
            cfl = METHOD_CFL[method]
            for backend in args.backends:
                for resolution in args.resolutions:
                    key = (case_name, method, backend, resolution)
                    if key in done:
                        print(
                            f"\n[report2-discontinuous] skipping completed "
                            f"case={case_name} method={method} "
                            f"backend={backend} n={resolution}"
                        )
                        continue
                    command = [
                        str(executable),
                        f"--backend={backend}",
                        f"--problem={case_name}",
                        f"--method={method}",
                        f"--nx={resolution}",
                        f"--ny={resolution}",
                        f"--boundary={case['boundary']}",
                        f"--final_time={final_times[case_name]}",
                        f"--cfl={cfl}",
                        "--output=none",
                    ]
                    environment = os.environ.copy()
                    if backend == "openmp":
                        environment["OMP_NUM_THREADS"] = str(args.openmp_threads)
                        environment["OMP_DYNAMIC"] = "FALSE"
                        environment.setdefault("OMP_PROC_BIND", "close")
                        environment.setdefault("OMP_PLACES", "cores")

                    print(
                        f"\n[report2-discontinuous] case={case_name} "
                        f"method={method} backend={backend} "
                        f"n={resolution} cfl={cfl}"
                    )
                    diagnostics_output = ""
                    for index in range(args.warmup):
                        elapsed, diagnostics_output = timed_run(command, environment)
                        print(f"  warmup {index + 1}: {elapsed:.6f} s")

                    samples: list[float] = []
                    for index in range(args.repeats):
                        elapsed, diagnostics_output = timed_run(command, environment)
                        samples.append(elapsed)
                        print(f"  repeat {index + 1}: {elapsed:.6f} s")

                    mean_seconds = statistics.fmean(samples)
                    stdev_seconds = (
                        statistics.stdev(samples) if len(samples) > 1 else 0.0
                    )
                    append_row(
                        output_path,
                        {
                            "problem": case_name,
                            "method": method,
                            "backend": backend,
                            "nx": resolution,
                            "ny": resolution,
                            "cfl": cfl,
                            "final_time": final_times[case_name],
                            "openmp_threads": (
                                args.openmp_threads if backend == "openmp" else ""
                            ),
                            "steps": parse_steps(diagnostics_output),
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
                    done.add(key)
                    print(f"  saved row to: {output_path}")

    print(f"\nsaved Report 2 discontinuous runtime CSV: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
