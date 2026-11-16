#!/usr/bin/env python3
"""
Run Stage B consistency checks for current high-order 1D smooth problems.

Current scope:
- problems: smooth_advection, gaussian_advection
- methods: weno2, weno3, ader3
- backends: serial, openmp

Rationale:
- ADER-3 is currently 1D, while WENO-3 also has a 2D path.
- The high-order 1D CUDA paths are not implemented.
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
SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
NX = 200
BUILD_DIR = "build-report2"


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    output_dir = root / "results" / "report2" / "backend_checks_1d"
    output_dir.mkdir(parents=True, exist_ok=True)
    executable = f"./{BUILD_DIR}/final_project"

    for case in CASES:
        for method in METHODS:
            serial_csv = (
                f"results/report2/backend_checks_1d/"
                f"{case['problem']}_serial_{method}.csv"
            )
            openmp_csv = (
                f"results/report2/backend_checks_1d/"
                f"{case['problem']}_openmp_{method}.csv"
            )

            base_args = [
                f"--problem={case['problem']}",
                f"--method={method}",
                f"--nx={NX}",
                "--boundary=periodic",
                f"--final_time={case['final_time']}",
                f"--cfl={SMOOTH_CFL[method]}",
            ]

            run(
                [executable, "--backend=serial", *base_args, f"--output={serial_csv}"],
                root,
            )
            run(
                [executable, "--backend=openmp", *base_args, f"--output={openmp_csv}"],
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
