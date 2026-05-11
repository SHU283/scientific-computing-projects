#!/usr/bin/env python3
"""
Generate missing representative Stage A artifacts for Report 1 methods:
- Godunov
- MUSCL
- WENO-2
"""

from __future__ import annotations

import subprocess
from pathlib import Path


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent.parent

    # 1D smooth representative runs
    for problem in ["smooth_advection", "gaussian_advection"]:
        for method in ["godunov", "muscl", "weno2"]:
            output_csv = f"results/stageA/{problem}_1d_{method}.csv"
            run(
                [
                    "./build/final_project",
                    "--backend=serial",
                    f"--problem={problem}",
                    f"--method={method}",
                    "--nx=400",
                    "--boundary=periodic",
                    "--final_time=0.3",
                    f"--output={output_csv}",
                ],
                root,
            )
            run(
                [
                    "python3",
                    str(root / "scripts" / "plot_field_1d.py"),
                    output_csv,
                    "--field",
                    "rho",
                    "--with-exact",
                ],
                root,
            )

    # 2D smooth representative MUSCL run, to pair with existing Godunov artifacts
    run(
        [
            "./build/final_project",
            "--backend=serial",
            "--problem=gaussian_advection_2d",
            "--method=muscl",
            "--nx=100",
            "--ny=100",
            "--final_time=0.1",
            "--output=results/stageA/gaussian_advection_2d_muscl.csv",
        ],
        root,
    )
    run(
        [
            "python3",
            str(root / "scripts" / "plot_field_2d.py"),
            "results/stageA/gaussian_advection_2d_muscl.csv",
            "--field=rho",
        ],
        root,
    )

    print("\nReport 1 Stage A representative artifacts finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
