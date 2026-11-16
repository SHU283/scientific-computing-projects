#!/usr/bin/env python3
"""Collect opt-in solver phase profiles for the smooth 2D HLLC methods."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("weno2", "weno3", "ader3")
BACKENDS = ("serial", "openmp", "cuda")
METHOD_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
PROFILE_PATTERN = re.compile(r"^PHASE_PROFILE\s+(?P<body>.+)$", re.MULTILINE)
PAIR_PATTERN = re.compile(r"(?P<key>[a-z_]+)=(?P<value>\S+)")
PROFILE_FIELDS = (
    "total_s",
    "setup_s",
    "cfl_s",
    "workspace_s",
    "flux_x_s",
    "flux_y_s",
    "spatial_fused_s",
    "rhs_assembly_s",
    "state_update_s",
    "admissibility_s",
    "transfer_s",
    "unattributed_s",
)
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
    "retries",
    "repeat",
    "wall_seconds",
    *PROFILE_FIELDS,
    "spatial_operator_s",
    "accounted_s",
)


def find_executable(build_dir: str) -> Path:
    candidate = ROOT / build_dir / "final_project"
    if candidate.exists():
        return candidate
    raise SystemExit(f"solver executable not found: {candidate}")


def parse_profile(output: str) -> dict[str, str]:
    match = PROFILE_PATTERN.search(output)
    if not match:
        raise RuntimeError(f"PHASE_PROFILE line was not found:\n{output}")
    values = {
        item.group("key"): item.group("value")
        for item in PAIR_PATTERN.finditer(match.group("body"))
    }
    missing = [field for field in (*PROFILE_FIELDS, "steps", "retries") if field not in values]
    if missing:
        raise RuntimeError(f"missing profile fields {missing}:\n{match.group(0)}")
    return values


def read_completed(path: Path) -> set[tuple[str, str, int, int]]:
    if not path.exists():
        return set()
    with path.open(newline="", encoding="utf-8") as stream:
        return {
            (row["method"], row["backend"], int(row["nx"]), int(row["repeat"]))
            for row in csv.DictReader(stream)
        }


def append_row(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(100, 600))
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--backends", nargs="+", choices=BACKENDS, default=BACKENDS)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--openmp-threads", type=int, default=16)
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--output",
        default="results/report2_hllc/performance/phase_profile/phase_profile_2d_hllc_lovelace_a30.csv",
    )
    args = parser.parse_args()

    executable = find_executable(args.build_dir)
    output_path = ROOT / args.output
    completed = read_completed(output_path) if args.resume else set()
    if output_path.exists() and not args.resume:
        output_path.unlink()

    for resolution in args.resolutions:
        for method in args.methods:
            for backend in args.backends:
                for repeat in range(1, args.repeats + 1):
                    key = (method, backend, resolution, repeat)
                    if key in completed:
                        print(f"[skip] method={method} backend={backend} n={resolution} repeat={repeat}")
                        continue

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
                    environment["FP_PHASE_PROFILE"] = "1"
                    if backend == "openmp":
                        environment["OMP_NUM_THREADS"] = str(args.openmp_threads)
                        environment["OMP_DYNAMIC"] = "FALSE"
                        environment.setdefault("OMP_PROC_BIND", "close")
                        environment.setdefault("OMP_PLACES", "cores")

                    print(
                        f"\n[phase-profile] method={method} backend={backend} "
                        f"n={resolution} repeat={repeat}"
                    )
                    start = time.perf_counter()
                    result = subprocess.run(
                        command,
                        cwd=ROOT,
                        text=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        env=environment,
                    )
                    wall_seconds = time.perf_counter() - start
                    if result.returncode != 0:
                        print(result.stdout)
                        raise subprocess.CalledProcessError(result.returncode, command)

                    profile = parse_profile(result.stdout)
                    numeric = {field: float(profile[field]) for field in PROFILE_FIELDS}
                    spatial = sum(
                        numeric[field]
                        for field in (
                            "flux_x_s",
                            "flux_y_s",
                            "spatial_fused_s",
                            "rhs_assembly_s",
                        )
                    )
                    accounted = numeric["total_s"] - numeric["unattributed_s"]
                    row: dict[str, object] = {
                        "problem": "gaussian_advection_2d",
                        "method": method,
                        "backend": backend,
                        "nx": resolution,
                        "ny": resolution,
                        "cfl": cfl,
                        "final_time": args.final_time,
                        "openmp_threads": args.openmp_threads if backend == "openmp" else "",
                        "steps": int(profile["steps"]),
                        "retries": int(profile["retries"]),
                        "repeat": repeat,
                        "wall_seconds": wall_seconds,
                        **numeric,
                        "spatial_operator_s": spatial,
                        "accounted_s": accounted,
                    }
                    append_row(output_path, row)
                    print(
                        f"  total={numeric['total_s']:.6f}s spatial={spatial:.6f}s "
                        f"cfl={numeric['cfl_s']:.6f}s update={numeric['state_update_s']:.6f}s"
                    )
                    print(f"  saved row: {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
