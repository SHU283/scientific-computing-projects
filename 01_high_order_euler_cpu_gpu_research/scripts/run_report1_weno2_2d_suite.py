#!/usr/bin/env python3
"""
Run the Report 1 WENO-2 2D smooth suite on gaussian_advection_2d.

Stages covered:
- Stage A: representative single run + convergence
- Stage B: serial/openmp consistency
- Stage C: runtime scaling for serial/openmp
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path


RESOLUTIONS = [50, 100, 200]
BACKENDS = ["serial", "openmp"]


def run(cmd: list[str], cwd: Path) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)
    print()


def timed_run(cmd: list[str], cwd: Path) -> float:
    start = time.perf_counter()
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return time.perf_counter() - start


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--runtime-output", default="results/stageC/runtime_scaling_2d_weno2.csv")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    exe = "./build/final_project"

    # Stage A representative single run.
    rep_csv = "results/stageA/gaussian_advection_2d_weno2.csv"
    run(
        [
            exe,
            "--backend=serial",
            "--problem=gaussian_advection_2d",
            "--method=weno2",
            "--nx=100",
            "--ny=100",
            "--boundary=periodic",
            "--final_time=0.1",
            f"--output={rep_csv}",
        ],
        root,
    )
    run(
        [
            "python3",
            "scripts/plot_field_2d.py",
            rep_csv,
            "--field=rho",
        ],
        root,
    )
    run(
        [
            exe,
            "--study=convergence",
            "--problem=gaussian_advection_2d",
            "--method=weno2",
        ],
        root,
    )

    # Stage B consistency.
    serial_csv = "results/stageB/gaussian_advection_2d_serial_weno2.csv"
    openmp_csv = "results/stageB/gaussian_advection_2d_openmp_weno2.csv"
    for backend, output_csv in [("serial", serial_csv), ("openmp", openmp_csv)]:
        run(
            [
                exe,
                f"--backend={backend}",
                "--problem=gaussian_advection_2d",
                "--method=weno2",
                "--nx=100",
                "--ny=100",
                "--boundary=periodic",
                "--final_time=0.1",
                f"--output={output_csv}",
            ],
            root,
        )
    run(["python3", "scripts/compare_csv_results.py", serial_csv, openmp_csv], root)

    # Stage C runtime scaling.
    rows: list[dict[str, str | int | float]] = []
    for n in RESOLUTIONS:
        for backend in BACKENDS:
            output_csv = f"results/stageC/runtime_scaling_2d_{backend}_weno2_n{n}.csv"
            cmd = [
                exe,
                f"--backend={backend}",
                "--problem=gaussian_advection_2d",
                "--method=weno2",
                f"--nx={n}",
                f"--ny={n}",
                "--final_time=0.1",
                f"--output={output_csv}",
            ]
            print(f"[weno2-2d-runtime] n={n} backend={backend}")
            for warmup in range(args.warmup):
                elapsed = timed_run(cmd, root)
                print(f"  warmup {warmup + 1}: {elapsed:.6f} s")
            samples: list[float] = []
            for repeat in range(args.repeats):
                elapsed = timed_run(cmd, root)
                samples.append(elapsed)
                print(f"  repeat {repeat + 1}: {elapsed:.6f} s")

            rows.append(
                {
                    "problem": "gaussian_advection_2d",
                    "backend": backend,
                    "method": "weno2",
                    "nx": n,
                    "ny": n,
                    "warmup": args.warmup,
                    "repeats": args.repeats,
                    "best_seconds": min(samples),
                    "mean_seconds": statistics.fmean(samples),
                    "median_seconds": statistics.median(samples),
                }
            )
            print()

    runtime_path = root / args.runtime_output
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    with runtime_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "problem",
                "backend",
                "method",
                "nx",
                "ny",
                "warmup",
                "repeats",
                "best_seconds",
                "mean_seconds",
                "median_seconds",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"saved WENO-2 2D runtime csv: {runtime_path}")
    print("WENO-2 2D suite finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
