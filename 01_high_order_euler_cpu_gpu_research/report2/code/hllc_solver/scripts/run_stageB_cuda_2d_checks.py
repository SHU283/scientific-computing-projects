#!/usr/bin/env python3
"""Run Stage B serial/CUDA consistency checks for smooth 2D cases.

The campaign covers Gaussian advection with Godunov and MUSCL. It exports both
backends and delegates field comparison to compare_csv_results.py.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

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

    for method in METHODS:
        serial_csv = f"results/stageB/gaussian_advection_2d_serial_{method}.csv"
        cuda_csv = f"results/stageB/gaussian_advection_2d_cuda_{method}.csv"

        run(
            [
                args.serial_exe,
                "--backend=serial",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                "--nx=50",
                "--ny=50",
                "--final_time=0.1",
                f"--output={serial_csv}",
            ],
            root,
        )

        run(
            [
                args.cuda_exe,
                "--backend=cuda",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                "--nx=50",
                "--ny=50",
                "--final_time=0.1",
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

    print("\nStage B 2D CUDA consistency check finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

