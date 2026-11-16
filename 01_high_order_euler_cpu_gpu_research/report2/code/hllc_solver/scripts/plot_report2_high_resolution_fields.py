#!/usr/bin/env python3
"""Render publication-ready Report 2 highest-resolution field comparisons."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from compare_fields_1d import exact_rho
from plot_sod_report1_points_exact import sod_exact_rho


METHODS = ("weno2", "weno3", "ader3")
LABELS = {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}
COLORS = {"weno2": "#d97706", "weno3": "#0f766e", "ader3": "#6d28d9"}
MARKERS = {"weno2": "s", "weno3": "^", "ader3": "D"}


def load_1d(path: Path) -> tuple[list[float], list[float]]:
    x: list[float] = []
    rho: list[float] = []
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            x.append(float(row["x"]))
            rho.append(float(row["rho"]))
    return x, rho


def plot_1d_case(
    input_dir: Path,
    output_dir: Path,
    problem: str,
    resolution: int,
    final_time: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["text.parse_math"] = False

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=300)
    for method in METHODS:
        path = input_dir / f"{problem}_{method}_n{resolution}.csv"
        if not path.exists():
            raise SystemExit(f"missing highest-resolution field: {path}")
        x, rho = load_1d(path)
        ax.plot(
            x,
            rho,
            color=COLORS[method],
            linewidth=2.1,
            marker=MARKERS[method],
            markersize=4.8,
            markevery=max(len(x) // 20, 1),
            label=LABELS[method],
        )

    exact_x = [index / 2000.0 for index in range(2001)]
    if problem == "sod":
        exact = [sod_exact_rho(value, final_time) for value in exact_x]
        title = f"Sod shock tube density at t = {final_time:g}, n = {resolution}"
    else:
        exact = [
            exact_rho(problem, value, final_time, 0.73, 0.0, 1.0)
            for value in exact_x
        ]
        display = (
            "Sinusoidal advection"
            if problem == "smooth_advection"
            else "Gaussian advection"
        )
        title = f"{display} density at t = {final_time:g}, n = {resolution}"

    ax.plot(exact_x, exact, "k--", linewidth=2.1, label="Exact")
    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("Density", fontsize=16)
    ax.set_title(title, fontsize=18, pad=12)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=13, ncol=2)
    fig.tight_layout()

    output = output_dir / f"{problem}_methods_n{resolution}.rho.png"
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


def plot_2d_cases(root: Path, input_dir: Path, output_dir: Path, resolution: int) -> None:
    cases = (
        ("gaussian_advection_2d", "2D Gaussian advection density at t = 0.1"),
        ("quadrant_riemann_2d", "2D quadrant Riemann density at t = 0.2"),
        ("shock_bubble_2d", "2D shock-bubble density at t = 0.2"),
    )
    comparison_script = root / "scripts" / "compare_fields_2d_multi.py"

    for problem, title in cases:
        available_methods = [
            method
            for method in METHODS
            if (input_dir / f"{problem}_{method}_n{resolution}.csv").exists()
        ]
        paths = [
            input_dir / f"{problem}_{method}_n{resolution}.csv"
            for method in available_methods
        ]
        if len(paths) < 2:
            print(
                f"skipping {problem}; fewer than two method fields are available",
                file=sys.stderr,
            )
            continue

        output = output_dir / f"{problem}_methods_n{resolution}.rho.png"
        subprocess.run(
            [
                sys.executable,
                str(comparison_script),
                *(str(path) for path in paths),
                "--field",
                "rho",
                "--suptitle",
                f"{title}, {resolution} x {resolution}",
                "--total-figsize",
                "16.8" if len(paths) == 3 else "12.4",
                "6.0",
                "--dpi",
                "300",
                "--title-fontsize",
                "22",
                "--label-fontsize",
                "18",
                "--tick-fontsize",
                "15",
                "--colorbar-fontsize",
                "18",
                "--output",
                str(output),
            ],
            cwd=root,
            check=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-dir", default="results/report2/high_resolution_fields"
    )
    parser.add_argument("--resolution-1d", type=int, default=800)
    parser.add_argument("--resolution-2d", type=int, default=600)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    base_dir = root / args.base_dir
    input_1d = base_dir / "1d"
    input_2d = base_dir / "2d"
    output_dir = base_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_1d_case(
        input_1d, output_dir, "smooth_advection", args.resolution_1d, 0.3
    )
    plot_1d_case(
        input_1d, output_dir, "gaussian_advection", args.resolution_1d, 0.3
    )
    plot_1d_case(input_1d, output_dir, "sod", args.resolution_1d, 0.2)
    plot_2d_cases(root, input_2d, output_dir, args.resolution_2d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
