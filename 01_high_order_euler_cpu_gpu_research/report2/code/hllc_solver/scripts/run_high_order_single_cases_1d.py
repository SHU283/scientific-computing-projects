#!/usr/bin/env python3
"""
Generate representative 1D single-case CSVs and figures for high-order methods.

This fills the Stage A / Stage D gap where some methods only had convergence
outputs but no dedicated representative single-run CSVs for plotting.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


PROBLEMS = ["smooth_advection", "gaussian_advection"]
METHODS = ["muscl", "weno2", "weno3", "ader3"]
REPORT2_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
NX = 400
FINAL_TIME = 0.3


def run(cmd: list[str], cwd: Path) -> None:
    print(f"\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> int:
    root = Path(__file__).resolve().parent.parent

    for problem in PROBLEMS:
        csv_paths: list[str] = []
        for method in METHODS:
            output_csv = f"results/stageA/{problem}_1d_{method}.csv"
            csv_paths.append(output_csv)
            command = [
                    "./build/final_project",
                    "--backend=serial",
                    f"--problem={problem}",
                    f"--method={method}",
                    f"--nx={NX}",
                    "--boundary=periodic",
                    f"--final_time={FINAL_TIME}",
                    f"--output={output_csv}",
                ]
            if method in REPORT2_CFL:
                command.insert(-1, f"--cfl={REPORT2_CFL[method]}")
            run(command, root)
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

        run(
            [
                "python3",
                str(root / "scripts" / "compare_fields_1d.py"),
                *csv_paths,
                "--field",
                "rho",
            ],
            root,
        )

    print("\nHigh-order 1D representative cases and figures finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
