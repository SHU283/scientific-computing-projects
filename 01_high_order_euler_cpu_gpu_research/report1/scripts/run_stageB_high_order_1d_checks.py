#!/usr/bin/env python3
"""
Run Stage B consistency checks for current high-order 1D smooth problems.

Current scope:
- problems: smooth_advection, gaussian_advection
- methods: weno2, weno3, ader3
- backends: serial, openmp

Rationale:
- These high-order implementations are currently 1D-only.
- CUDA paths are not implemented for them yet.
- OpenMP is still worth checking for result consistency, even if the
  current implementation is mostly shared code.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


CASES = [
    {"problem": "smooth_advection", "final_time": 0.3},
    {"problem": "gaussian_advection", "final_time": 0.3},
]

METHODS = ["weno2", "weno3", "ader3"]
NX = 200


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent.parent

    for case in CASES:
        for method in METHODS:
            serial_csv = f"results/stageB/{case['problem']}_serial_{method}.csv"
            openmp_csv = f"results/stageB/{case['problem']}_openmp_{method}.csv"

            base_args = [
                f"--problem={case['problem']}",
                f"--method={method}",
                f"--nx={NX}",
                "--boundary=periodic",
                f"--final_time={case['final_time']}",
            ]

            run(
                ["./build/final_project", "--backend=serial", *base_args, f"--output={serial_csv}"],
                root,
            )
            run(
                ["./build/final_project", "--backend=openmp", *base_args, f"--output={openmp_csv}"],
                root,
            )
            run(
                ["python3", str(root / "scripts" / "compare_csv_results.py"), serial_csv, openmp_csv],
                root,
            )

    print("\nStage B high-order 1D consistency checks finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

