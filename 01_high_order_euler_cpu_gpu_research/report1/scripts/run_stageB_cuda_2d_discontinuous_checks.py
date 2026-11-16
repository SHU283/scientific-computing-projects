#!/usr/bin/env python3
"""Run Stage B CUDA consistency checks for discontinuous 2D cases.

The selected quadrant Riemann and shock-bubble runs exercise the GPU backend
on sharper flow features before comparing the output against serial results.
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
    parser.add_argument("--serial-exe", default="./build/final_project")
    parser.add_argument("--cuda-exe", default="./build-cuda/final_project")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    compare_script = root / "scripts" / "compare_csv_results.py"

    for problem, nx, ny, final_time in CASES:
        for method in METHODS:
            serial_csv = f"results/stageB/{problem}_serial_{method}.csv"
            cuda_csv = f"results/stageB/{problem}_cuda_{method}.csv"

            run(
                [
                    args.serial_exe,
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
                    args.cuda_exe,
                    "--backend=cuda",
                    f"--problem={problem}",
                    f"--method={method}",
                    f"--nx={nx}",
                    f"--ny={ny}",
                    f"--final_time={final_time}",
                    f"--output={cuda_csv}",
                ],
                root,
            )

            run(
                [
                    args.python,
                    str(compare_script),
                    serial_csv,
                    cuda_csv,
                ],
                root,
            )

    print("\nStage B 2D discontinuous CUDA consistency check finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


