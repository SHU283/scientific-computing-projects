#!/usr/bin/env python3
"""Run small serial/OpenMP cases and compare their exported 2D fields."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("weno2", "weno3", "ader3")
CASES = {
    "smooth": {
        "problem": "gaussian_advection_2d",
        "n": 50,
        "boundary": "periodic",
        "final_time": 0.02,
    },
    "quadrant": {
        "problem": "quadrant_riemann_2d",
        "n": 30,
        "boundary": "outflow",
        "final_time": 0.01,
    },
    "shock": {
        "problem": "shock_bubble_2d",
        "n": 30,
        "boundary": "outflow",
        "final_time": 0.005,
    },
}
SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
DISCONTINUOUS_CFL = {"weno2": 0.2, "weno3": 0.2, "ader3": 0.14}
FIELDS = ("rho", "u", "v", "p", "energy")
STEP_PATTERN = re.compile(r"time steps\s+:\s+(\d+)")


def run_case(
    executable: Path,
    output_dir: Path,
    case_name: str,
    method: str,
    backend: str,
    openmp_threads: int,
) -> tuple[Path, int]:
    case = CASES[case_name]
    cfl = SMOOTH_CFL[method] if case_name == "smooth" else DISCONTINUOUS_CFL[method]
    output_path = output_dir / f"{case_name}_{method}_{backend}.csv"
    command = [
        str(executable),
        f"--backend={backend}",
        f"--problem={case['problem']}",
        f"--method={method}",
        f"--nx={case['n']}",
        f"--ny={case['n']}",
        f"--boundary={case['boundary']}",
        f"--final_time={case['final_time']}",
        f"--cfl={cfl}",
        f"--output={output_path}",
    ]
    environment = os.environ.copy()
    if backend == "openmp":
        environment.update(
            {
                "OMP_NUM_THREADS": str(openmp_threads),
                "OMP_DYNAMIC": "FALSE",
                "OMP_PROC_BIND": "close",
                "OMP_PLACES": "cores",
            }
        )
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    match = STEP_PATTERN.search(completed.stdout)
    if not match:
        raise RuntimeError(f"time-step count missing from output:\n{completed.stdout}")
    return output_path, int(match.group(1))


def compare_fields(serial_path: Path, openmp_path: Path) -> dict[str, float]:
    maxima = {field: 0.0 for field in FIELDS}
    with serial_path.open(newline="", encoding="utf-8") as serial_stream, openmp_path.open(
        newline="", encoding="utf-8"
    ) as openmp_stream:
        serial_reader = csv.DictReader(serial_stream)
        openmp_reader = csv.DictReader(openmp_stream)
        serial_rows = list(serial_reader)
        openmp_rows = list(openmp_reader)
    if len(serial_rows) != len(openmp_rows):
        raise RuntimeError(
            f"row-count mismatch: {serial_path}={len(serial_rows)}, "
            f"{openmp_path}={len(openmp_rows)}"
        )
    for serial_row, openmp_row in zip(serial_rows, openmp_rows, strict=True):
        if (serial_row["i"], serial_row["j"]) != (openmp_row["i"], openmp_row["j"]):
            raise RuntimeError("cell ordering differs between serial and OpenMP outputs")
        for field in FIELDS:
            maxima[field] = max(
                maxima[field],
                abs(float(serial_row[field]) - float(openmp_row[field])),
            )
    return maxima


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--executable",
        type=Path,
        default=ROOT / "build-hllc-local" / "final_project",
    )
    parser.add_argument("--openmp-threads", type=int, default=4)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results" / "_alignment_audit",
    )
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    executable = args.executable.resolve()
    if not executable.exists():
        raise SystemExit(f"executable not found: {executable}")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    failed = False
    for case_name in CASES:
        for method in METHODS:
            serial_path, serial_steps = run_case(
                executable, output_dir, case_name, method, "serial", args.openmp_threads
            )
            openmp_path, openmp_steps = run_case(
                executable, output_dir, case_name, method, "openmp", args.openmp_threads
            )
            maxima = compare_fields(serial_path, openmp_path)
            max_delta = max(maxima.values())
            status = "PASS"
            if serial_steps != openmp_steps or max_delta > args.tolerance:
                status = "FAIL"
                failed = True
            print(
                f"{status} case={case_name} method={method} "
                f"steps={serial_steps}/{openmp_steps} max_delta={max_delta:.3e} "
                f"fields={maxima}"
            )

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
