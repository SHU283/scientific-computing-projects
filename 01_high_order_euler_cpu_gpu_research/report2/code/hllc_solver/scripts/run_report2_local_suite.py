#!/usr/bin/env python3
"""Run the reproducible local Report 2 correctness and convergence suite."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(command: list[str], root: Path) -> None:
    print(f"\n$ {' '.join(command)}", flush=True)
    subprocess.run(command, cwd=root, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2")
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Use the existing executable without rebuilding it.",
    )
    parser.add_argument(
        "--full-2d",
        action="store_true",
        help="Extend the 2D convergence run from n=100 to n=200 and n=400.",
    )
    parser.add_argument(
        "--skip-1d",
        action="store_true",
        help="Skip the established 1D convergence campaign.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    python = sys.executable
    build_dir = root / args.build_dir
    executable = build_dir / "final_project"

    if not args.skip_build:
        run(["cmake", "--build", args.build_dir, "-j"], root)
    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    run(
        [
            python,
            "scripts/run_report2_local_smoke.py",
            "--build-dir",
            args.build_dir,
        ],
        root,
    )

    if not args.skip_1d:
        run(
            [
                python,
                "scripts/run_report2_validation_1d.py",
                "--build-dir",
                args.build_dir,
            ],
            root,
        )

    resolutions = ["25", "50", "100"]
    if args.full_2d:
        resolutions.extend(["200", "400"])
    unified_dir = "results/report2/validation_2d_unified"
    run(
        [
            python,
            "scripts/run_report2_validation_2d.py",
            "--build-dir",
            args.build_dir,
            "--output-dir",
            unified_dir,
            "--methods",
            "weno2",
            "weno3",
            "ader3",
            "--resolutions",
            *resolutions,
        ],
        root,
    )
    run(
        [
            python,
            "scripts/plot_report2_convergence_2d.py",
            f"{unified_dir}/report2_high_order_convergence_2d.csv",
            "--output",
            f"{unified_dir}/report2_high_order_convergence_2d.l1_rho.png",
        ],
        root,
    )

    print("\nReport 2 local suite: PASS")
    print(f"Unified 2D evidence: {root / unified_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
