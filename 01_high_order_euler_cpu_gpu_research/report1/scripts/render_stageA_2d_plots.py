#!/usr/bin/env python3
"""
Batch-render the core Stage A 2D figures with unified size and DPI.

This script covers the representative smooth-case density plots and a more
complete set of discontinuous-case plots so the results folder stays organized.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


DEFAULT_PLOTS = [
    # 2D smooth validation: density is the most informative field here.
    ("results/stageA/gaussian_advection_2d_godunov.csv", "rho"),
    ("results/stageA/gaussian_advection_2d_muscl.csv", "rho"),
    # 2D discontinuous validation: keep rho / p / energy for both methods.
    ("results/stageA/quadrant_riemann_2d_godunov.csv", "rho"),
    ("results/stageA/quadrant_riemann_2d_godunov.csv", "p"),
    ("results/stageA/quadrant_riemann_2d_godunov.csv", "energy"),
    ("results/stageA/quadrant_riemann_2d_muscl.csv", "rho"),
    ("results/stageA/quadrant_riemann_2d_muscl.csv", "p"),
    ("results/stageA/quadrant_riemann_2d_muscl.csv", "energy"),
    ("results/stageA/shock_bubble_2d_godunov.csv", "rho"),
    ("results/stageA/shock_bubble_2d_godunov.csv", "p"),
    ("results/stageA/shock_bubble_2d_godunov.csv", "energy"),
    ("results/stageA/shock_bubble_2d_muscl.csv", "rho"),
    ("results/stageA/shock_bubble_2d_muscl.csv", "p"),
    ("results/stageA/shock_bubble_2d_muscl.csv", "energy"),
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(9.0, 7.0))
    parser.add_argument("--interpolation", default="nearest")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    plot_script = root / "scripts" / "plot_field_2d.py"

    for csv_rel, field in DEFAULT_PLOTS:
        csv_path = root / csv_rel
        if not csv_path.exists():
            print(f"skip missing file: {csv_path}")
            continue

        cmd = [
            sys.executable,
            str(plot_script),
            str(csv_path),
            "--field",
            field,
            "--dpi",
            str(args.dpi),
            "--figsize",
            str(args.figsize[0]),
            str(args.figsize[1]),
            "--interpolation",
            args.interpolation,
        ]
        subprocess.run(cmd, check=True)

    print("all available Stage A 2D plots rendered.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


