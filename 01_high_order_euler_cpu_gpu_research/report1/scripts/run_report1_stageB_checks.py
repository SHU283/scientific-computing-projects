#!/usr/bin/env python3
"""
Fill the currently missing Report 1 Stage B consistency archives.

Focus:
- low-order methods: Godunov, MUSCL
- report-1 methods only
- missing serial/openmp consistency sets for 1D smooth, 1D discontinuous, and 2D smooth
"""

from __future__ import annotations

import subprocess
from pathlib import Path


ONE_D_CASES = [
    {"problem": "smooth_advection", "method": "godunov", "nx": 200, "final_time": 0.3, "boundary": "periodic"},
    {"problem": "smooth_advection", "method": "muscl", "nx": 200, "final_time": 0.3, "boundary": "periodic"},
    {"problem": "gaussian_advection", "method": "godunov", "nx": 200, "final_time": 0.3, "boundary": "periodic"},
    {"problem": "gaussian_advection", "method": "muscl", "nx": 200, "final_time": 0.3, "boundary": "periodic"},
    {"problem": "sod", "method": "godunov", "nx": 200, "final_time": 0.2, "boundary": "outflow"},
    {"problem": "sod", "method": "muscl", "nx": 200, "final_time": 0.2, "boundary": "outflow"},
]

TWO_D_CASES = [
    {"problem": "gaussian_advection_2d", "method": "godunov", "nx": 100, "ny": 100, "final_time": 0.1},
    {"problem": "gaussian_advection_2d", "method": "muscl", "nx": 100, "ny": 100, "final_time": 0.1},
]


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    compare_script = str(root / "scripts" / "compare_csv_results.py")

    for case in ONE_D_CASES:
        serial_csv = f"results/stageB/{case['problem']}_serial_{case['method']}.csv"
        openmp_csv = f"results/stageB/{case['problem']}_openmp_{case['method']}.csv"

        base_args = [
            f"--problem={case['problem']}",
            f"--method={case['method']}",
            f"--nx={case['nx']}",
            f"--boundary={case['boundary']}",
            f"--final_time={case['final_time']}",
        ]

        run(["./build/final_project", "--backend=serial", *base_args, f"--output={serial_csv}"], root)
        run(["./build/final_project", "--backend=openmp", *base_args, f"--output={openmp_csv}"], root)
        run(["python3", compare_script, serial_csv, openmp_csv], root)

    for case in TWO_D_CASES:
        serial_csv = f"results/stageB/{case['problem']}_serial_{case['method']}.csv"
        openmp_csv = f"results/stageB/{case['problem']}_openmp_{case['method']}.csv"

        base_args = [
            f"--problem={case['problem']}",
            f"--method={case['method']}",
            f"--nx={case['nx']}",
            f"--ny={case['ny']}",
            "--boundary=periodic",
            f"--final_time={case['final_time']}",
        ]

        run(["./build/final_project", "--backend=serial", *base_args, f"--output={serial_csv}"], root)
        run(["./build/final_project", "--backend=openmp", *base_args, f"--output={openmp_csv}"], root)
        run(["python3", compare_script, serial_csv, openmp_csv], root)

    print("\nReport 1 Stage B missing consistency checks finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

