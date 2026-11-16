#!/usr/bin/env python3
"""Benchmark the Report 2 methods on the smooth 2D Gaussian problem."""

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
DEFAULT_RESOLUTIONS = (25, 50, 100, 200, 600)
METHOD_CFL = {
    "weno2": 0.4,
    "weno3": 0.4,
    "ader3": 0.25,
}
ERROR_PATTERN = re.compile(
    r"2D density error L1\s+:\s+(?P<l1>\S+).*?"
    r"2D density error L2\s+:\s+(?P<l2>\S+).*?"
    r"2D density error Linf\s+:\s+(?P<linf>\S+)",
    re.DOTALL,
)
STEP_PATTERN = re.compile(r"time steps\s+:\s+(?P<steps>\d+)")


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
    environment: dict[str, str] | None,
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


def parse_diagnostics(output: str) -> tuple[int, float, float, float]:
    error_match = ERROR_PATTERN.search(output)
    step_match = STEP_PATTERN.search(output)
    if not error_match or not step_match:
        raise RuntimeError(f"could not parse solver diagnostics:\n{output}")
    return (
        int(step_match.group("steps")),
        float(error_match.group("l1")),
        float(error_match.group("l2")),
        float(error_match.group("linf")),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default=None)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--resolutions",
        type=int,
        nargs="+",
        default=DEFAULT_RESOLUTIONS,
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=METHODS,
    )
    parser.add_argument(
        "--backends",
        nargs="+",
        choices=BACKENDS,
        default=BACKENDS,
    )
    parser.add_argument("--openmp-threads", type=int, default=16)
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument(
        "--output",
        default=(
            "results/report2/performance/"
            "runtime_accuracy_2d_report2_lovelace_a30.csv"
        ),
    )
    args = parser.parse_args()

    executable = find_executable(args.build_dir)
    rows: list[dict[str, str | int | float]] = []

    for method in args.methods:
        for backend in args.backends:
            for resolution in args.resolutions:
                cfl = METHOD_CFL[method]
                command = [
                    str(executable),
                    f"--backend={backend}",
                    "--problem=gaussian_advection_2d",
                    f"--method={method}",
                    f"--nx={resolution}",
                    f"--ny={resolution}",
                    "--boundary=periodic",
                    f"--final_time={args.final_time}",
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
                    f"\n[report2-runtime] method={method} "
                    f"backend={backend} n={resolution} cfl={cfl}"
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

                steps, l1, l2, linf = parse_diagnostics(diagnostics_output)
                rows.append(
                    {
                        "problem": "gaussian_advection_2d",
                        "method": method,
                        "backend": backend,
                        "nx": resolution,
                        "ny": resolution,
                        "cfl": cfl,
                        "final_time": args.final_time,
                        "openmp_threads": (
                            args.openmp_threads if backend == "openmp" else ""
                        ),
                        "steps": steps,
                        "warmup": args.warmup,
                        "repeats": args.repeats,
                        "best_seconds": min(samples),
                        "mean_seconds": statistics.fmean(samples),
                        "median_seconds": statistics.median(samples),
                        "l1_rho": l1,
                        "l2_rho": l2,
                        "linf_rho": linf,
                    }
                )

    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nsaved Report 2 runtime/accuracy CSV: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
