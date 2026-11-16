#!/usr/bin/env python3
"""Generate the highest-resolution 1D fields used in Report 2 figures."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
CASES = (
    ("smooth_advection", "periodic", 0.3, "smooth"),
    ("gaussian_advection", "periodic", 0.3, "smooth"),
    ("sod", "outflow", 0.2, "discontinuous"),
)
SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
DISCONTINUOUS_CFL = {"weno2": 0.2, "weno3": 0.2, "ader3": 0.14}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2")
    parser.add_argument(
        "--output-dir", default="results/report2/high_resolution_fields/1d"
    )
    parser.add_argument("--resolution", type=int, default=800)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    executable = root / args.build_dir / "final_project"
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    for problem, boundary, final_time, case_type in CASES:
        for method in METHODS:
            cfl = (
                SMOOTH_CFL[method]
                if case_type == "smooth"
                else DISCONTINUOUS_CFL[method]
            )
            output = output_dir / f"{problem}_{method}_n{args.resolution}.csv"
            command = [
                str(executable),
                "--backend=serial",
                f"--problem={problem}",
                f"--method={method}",
                f"--nx={args.resolution}",
                f"--boundary={boundary}",
                f"--final_time={final_time}",
                f"--cfl={cfl}",
                f"--output={output.relative_to(root).as_posix()}",
            ]
            print(
                f"[report2-high-resolution-1d] problem={problem} "
                f"method={method} n={args.resolution} cfl={cfl}"
            )
            subprocess.run(command, cwd=root, check=True)

    print(f"saved highest-resolution 1D fields below: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
