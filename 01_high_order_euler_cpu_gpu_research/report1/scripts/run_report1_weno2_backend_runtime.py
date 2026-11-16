#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

CASES = [
    {
        "problem": "smooth_advection",
        "nx": 200,
        "boundary": "periodic",
        "final_time": 0.3,
    },
    {
        "problem": "gaussian_advection",
        "nx": 200,
        "boundary": "periodic",
        "final_time": 0.3,
    },
    {
        "problem": "sod",
        "nx": 200,
        "boundary": "outflow",
        "final_time": 0.2,
    },
    {
        "problem": "gaussian_advection_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "periodic",
        "final_time": 0.1,
    },
    {
        "problem": "quadrant_riemann_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "outflow",
        "final_time": 0.2,
    },
    {
        "problem": "shock_bubble_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "outflow",
        "final_time": 0.05,
    },
]

BACKENDS = ["serial", "openmp", "cuda"]


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


def solver_cmd(executable: Path, backend: str, case: dict[str, object], output: str) -> list[str]:
    cmd = [
        str(executable),
        f"--backend={backend}",
        f"--problem={case['problem']}",
        "--method=weno2",
        f"--nx={case['nx']}",
        f"--boundary={case['boundary']}",
        f"--final_time={case['final_time']}",
        f"--output={output}",
    ]
    if "ny" in case:
        cmd.insert(5, f"--ny={case['ny']}")
    return cmd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output", default="results/stageC/runtime_weno2_backend_alignment.csv")
    args = parser.parse_args()

    executable = find_executable()
    if "build-cuda" not in executable.as_posix():
        raise SystemExit("Please use a CUDA build first: cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON")

    output_path = ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int | float]] = []
    for case in CASES:
        for backend in BACKENDS:
            output_csv = f"results/stageC/bench_{case['problem']}_{backend}_weno2.csv"
            cmd = solver_cmd(executable, backend, case, output_csv)

            print(f"\n[weno2-backend-runtime] problem={case['problem']} backend={backend}")
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
                    "problem": case["problem"],
                    "backend": backend,
                    "method": "weno2",
                    "nx": int(case["nx"]),
                    "ny": int(case.get("ny", 1)),
                    "boundary": str(case["boundary"]),
                    "final_time": float(case["final_time"]),
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
                "boundary",
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

    print(f"\nsaved WENO-2 backend runtime csv: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

