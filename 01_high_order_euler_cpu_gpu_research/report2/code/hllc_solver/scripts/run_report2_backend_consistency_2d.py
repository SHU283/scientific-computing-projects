#!/usr/bin/env python3
"""Run representative Report 2 2D backend-consistency checks."""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("weno2", "weno3", "ader3")
BACKENDS = ("serial", "openmp", "cuda")
FIELDS = ("rho", "u", "v", "p", "energy")
CASES = {
    "gaussian_advection_2d": {
        "boundary": "periodic",
        "final_time": 0.1,
        "smooth": True,
    },
    "quadrant_riemann_2d": {
        "boundary": "outflow",
        "final_time": 0.2,
        "smooth": False,
    },
    "shock_bubble_2d": {
        "boundary": "outflow",
        "final_time": 0.2,
        "smooth": False,
    },
}
SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
DISCONTINUOUS_CFL = {"weno2": 0.2, "weno3": 0.2, "ader3": 0.14}


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
    raise SystemExit("No final_project executable found.")


def run_solver(
    executable: Path,
    *,
    case_name: str,
    method: str,
    backend: str,
    resolution: int,
    output_path: Path,
    openmp_threads: int,
) -> None:
    case = CASES[case_name]
    cfl = SMOOTH_CFL[method] if case["smooth"] else DISCONTINUOUS_CFL[method]
    command = [
        str(executable),
        f"--backend={backend}",
        f"--problem={case_name}",
        f"--method={method}",
        f"--nx={resolution}",
        f"--ny={resolution}",
        f"--boundary={case['boundary']}",
        f"--final_time={case['final_time']}",
        f"--cfl={cfl}",
        f"--output={output_path}",
    ]
    environment = os.environ.copy()
    if backend == "openmp":
        environment["OMP_NUM_THREADS"] = str(openmp_threads)
        environment["OMP_DYNAMIC"] = "FALSE"
        environment.setdefault("OMP_PROC_BIND", "close")
        environment.setdefault("OMP_PLACES", "cores")

    print(
        f"[backend-consistency] case={case_name} method={method} "
        f"backend={backend} n={resolution}"
    )
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
    )
    if completed.returncode != 0:
        print(completed.stdout)
        raise subprocess.CalledProcessError(
            completed.returncode,
            command,
            output=completed.stdout,
        )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def min_field(rows: list[dict[str, str]], field: str) -> float:
    return min(float(row[field]) for row in rows)


def compare_rows(
    reference_rows: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
) -> dict[str, float]:
    if len(reference_rows) != len(candidate_rows):
        raise RuntimeError(
            f"row count mismatch: {len(reference_rows)} vs {len(candidate_rows)}"
        )
    differences: dict[str, float] = {}
    for field in FIELDS:
        max_difference = max(
            abs(float(reference[field]) - float(candidate[field]))
            for reference, candidate in zip(reference_rows, candidate_rows)
        )
        differences[field] = max_difference
    return differences


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default=None)
    parser.add_argument("--resolution", type=int, default=100)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--cases", nargs="+", choices=tuple(CASES), default=tuple(CASES))
    parser.add_argument("--openmp-threads", type=int, default=16)
    parser.add_argument(
        "--output-dir",
        default="results/report2/backend_consistency_2d_lovelace",
    )
    args = parser.parse_args()

    executable = find_executable(args.build_dir)
    output_dir = ROOT / args.output_dir
    fields_dir = output_dir / "fields"
    fields_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str | int | float]] = []
    for case_name in args.cases:
        for method in args.methods:
            case = CASES[case_name]
            cfl = SMOOTH_CFL[method] if case["smooth"] else DISCONTINUOUS_CFL[method]
            field_paths: dict[str, Path] = {}
            for backend in BACKENDS:
                output_path = (
                    fields_dir
                    / f"{case_name}_{method}_{backend}_n{args.resolution}.csv"
                )
                run_solver(
                    executable,
                    case_name=case_name,
                    method=method,
                    backend=backend,
                    resolution=args.resolution,
                    output_path=output_path,
                    openmp_threads=args.openmp_threads,
                )
                field_paths[backend] = output_path

            serial_rows = read_rows(field_paths["serial"])
            for backend in BACKENDS:
                rows = read_rows(field_paths[backend])
                if backend == "serial":
                    differences = {field: 0.0 for field in FIELDS}
                else:
                    differences = compare_rows(serial_rows, rows)
                finite = all(
                    math.isfinite(float(row[field]))
                    for row in rows
                    for field in FIELDS
                )
                summary_rows.append(
                    {
                        "problem": case_name,
                        "method": method,
                        "backend": backend,
                        "reference_backend": "serial",
                        "nx": args.resolution,
                        "ny": args.resolution,
                        "cfl": cfl,
                        "final_time": case["final_time"],
                        "openmp_threads": (
                            args.openmp_threads if backend == "openmp" else ""
                        ),
                        "rows": len(rows),
                        "finite": int(finite),
                        "min_rho": min_field(rows, "rho"),
                        "min_p": min_field(rows, "p"),
                        **{
                            f"max_abs_{field}_vs_serial": differences[field]
                            for field in FIELDS
                        },
                    }
                )

    summary_path = output_dir / "backend_consistency_2d_report2_lovelace_a30.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"saved backend consistency summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
