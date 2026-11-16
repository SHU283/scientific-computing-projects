#!/usr/bin/env python3
"""Run Stage B serial/OpenMP checks for discontinuous 2D problems.

The campaign covers the quadrant Riemann and shock-bubble cases, then compares
the exported serial and OpenMP fields with compare_csv_results.py.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


CASES = [
    ("quadrant_riemann_2d", 100, 100, 0.2),
    ("shock_bubble_2d", 100, 100, 0.05),
]

METHODS = ["godunov", "muscl"]


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", default="./build/final_project")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    compare_script = root / "scripts" / "compare_csv_results.py"

    for problem, nx, ny, final_time in CASES:
        for method in METHODS:
            serial_csv = f"results/stageB/{problem}_serial_{method}.csv"
            openmp_csv = f"results/stageB/{problem}_openmp_{method}.csv"

            run(
                [
                    args.exe,
                    "--backend=serial",
                    f"--problem={problem}",
                    f"--method={method}",
                    f"--nx={nx}",
                    f"--ny={ny}",
                    f"--final_time={final_time}",
                    f"--output={serial_csv}",
                ],
                root,
            )

            run(
                [
                    args.exe,
                    "--backend=openmp",
                    f"--problem={problem}",
                    f"--method={method}",
                    f"--nx={nx}",
                    f"--ny={ny}",
                    f"--final_time={final_time}",
                    f"--output={openmp_csv}",
                ],
                root,
            )

            run(
                [
                    args.python,
                    str(compare_script),
                    serial_csv,
                    openmp_csv,
                ],
                root,
            )

    print("\nStage B discontinuous consistency checks finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

