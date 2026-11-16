#!/usr/bin/env python3
"""Incremental CUDA-only high-resolution smooth sweep for Report 2.

This is intentionally separate from the canonical all-backend timing script.
It writes one CSV row immediately after each completed run so that long
high-resolution campaigns keep partial results if a later case fails.
"""

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
DEFAULT_RESOLUTIONS = (2400, 3000, 3600, 4200, 4600, 5000)
METHOD_CFL = {
    "weno2": 0.4,
    "weno3": 0.4,
    "ader3": 0.25,
}
FIELDNAMES = (
    "problem",
    "method",
    "backend",
    "nx",
    "ny",
    "cfl",
    "final_time",
    "steps",
    "warmup",
    "repeats",
    "best_seconds",
    "mean_seconds",
    "median_seconds",
    "l1_rho",
    "l2_rho",
    "linf_rho",
)
FAILURE_FIELDNAMES = (
    "problem",
    "method",
    "backend",
    "nx",
    "ny",
    "cfl",
    "final_time",
    "returncode",
    "message_tail",
)
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
    raise SystemExit("No CUDA-enabled final_project executable found.")


def timed_run(command: list[str], environment: dict[str, str]) -> tuple[float, str]:
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


def existing_keys(path: Path) -> set[tuple[str, int]]:
    if not path.exists():
        return set()
    with path.open(newline="", encoding="utf-8") as stream:
        return {
            (row["method"], int(row["nx"]))
            for row in csv.DictReader(stream)
        }


def append_row(path: Path, row: dict[str, str | int | float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()


def append_failure(path: Path, row: dict[str, str | int | float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FAILURE_FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default=None)
    parser.add_argument("--warmup", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=1)
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
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument(
        "--output",
        default=(
            "results/report2/performance/gpu_limit/"
            "runtime_accuracy_2d_report2_cuda_gpu_limit_extension_lovelace_a30.csv"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument(
        "--failure-output",
        default=(		
            "results/report2_hllc/performance/gpu_limit/"
            "runtime_accuracy_2d_report2_hllc_cuda_gpu_limit_extension_failures_lovelace_a30.csv"
        ),
    )
    args = parser.parse_args()

    executable = find_executable(args.build_dir)
    output_path = ROOT / args.output
    failure_path = ROOT / args.failure_output
    completed = existing_keys(output_path) if args.resume else set()

    environment = os.environ.copy()
    for method in args.methods:
        for resolution in args.resolutions:
            if (method, resolution) in completed:
                print(f"skip existing method={method} n={resolution}")
                continue
            cfl = METHOD_CFL[method]
            command = [
                str(executable),
                "--backend=cuda",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={resolution}",
                f"--ny={resolution}",
                "--boundary=periodic",
                f"--final_time={args.final_time}",
                f"--cfl={cfl}",
                "--output=none",
            ]
            print(f"\n[gpu-limit-extension] method={method} n={resolution} cfl={cfl}")

            try:
                diagnostics_output = ""
                for index in range(args.warmup):
                    elapsed, diagnostics_output = timed_run(command, environment)
                    print(f"  warmup {index + 1}: {elapsed:.6f} s")

                samples: list[float] = []
                for index in range(args.repeats):
                    elapsed, diagnostics_output = timed_run(command, environment)
                    samples.append(elapsed)
                    print(f"  repeat {index + 1}: {elapsed:.6f} s")
            except subprocess.CalledProcessError as exc:
                tail = "\n".join(str(exc.output or "").splitlines()[-12:])
                append_failure(
                    failure_path,
                    {
                        "problem": "gaussian_advection_2d",
                        "method": method,
                        "backend": "cuda",
                        "nx": resolution,
                        "ny": resolution,
                        "cfl": cfl,
                        "final_time": args.final_time,
                        "returncode": exc.returncode,
                        "message_tail": tail,
                    },
                )
                print(f"  failed; saved failure row to: {failure_path}")
                if args.keep_going:
                    continue
                raise

            steps, l1, l2, linf = parse_diagnostics(diagnostics_output)
            row = {
                "problem": "gaussian_advection_2d",
                "method": method,
                "backend": "cuda",
                "nx": resolution,
                "ny": resolution,
                "cfl": cfl,
                "final_time": args.final_time,
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
            append_row(output_path, row)
            print(f"  saved row to: {output_path}")

    print(f"\nsaved CUDA-only extension CSV: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
