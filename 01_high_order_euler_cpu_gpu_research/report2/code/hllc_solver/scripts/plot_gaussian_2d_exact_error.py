#!/usr/bin/env python3
"""
Plot absolute density error for 2D Gaussian advection CSV outputs.

The exact reference matches the manufactured 2D Gaussian advection problem in
src/problems/gaussian_advection_2d.cpp: periodic transport of a Gaussian density
perturbation with constant velocity components (adv_vel, 0.37).
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_rho(path: Path) -> tuple[list[list[float]], int, int]:
    rows: list[tuple[int, int, float]] = []
    max_i = -1
    max_j = -1
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            i = int(row["i"])
            j = int(row["j"])
            rows.append((i, j, float(row["rho"])))
            max_i = max(max_i, i)
            max_j = max(max_j, j)

    nx = max_i + 1
    ny = max_j + 1
    rho = [[0.0 for _ in range(nx)] for _ in range(ny)]
    for i, j, value in rows:
        rho[j][i] = value
    return rho, nx, ny


def wrap_periodic(x: float, x_min: float, x_max: float) -> float:
    length = x_max - x_min
    while x < x_min:
        x += length
    while x >= x_max:
        x -= length
    return x


def shortest_periodic_delta(x: float, center: float, length: float) -> float:
    delta = x - center
    if delta > 0.5 * length:
        delta -= length
    elif delta < -0.5 * length:
        delta += length
    return delta


def exact_rho_2d(
    i: int,
    j: int,
    nx: int,
    ny: int,
    final_time: float,
    adv_vel: float,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> float:
    length_x = x_max - x_min
    length_y = y_max - y_min
    dx = length_x / nx
    dy = length_y / ny
    x = x_min + (i + 0.5) * dx
    y = y_min + (j + 0.5) * dy
    x0 = wrap_periodic(x_min + 0.31 * length_x + adv_vel * final_time, x_min, x_max)
    y0 = wrap_periodic(y_min + 0.43 * length_y + 0.37 * final_time, y_min, y_max)
    sigma = 0.08 * length_x
    rx = shortest_periodic_delta(x, x0, length_x) / sigma
    ry = shortest_periodic_delta(y, y0, length_y) / sigma
    return 1.0 + 0.2 * math.exp(-0.5 * (rx * rx + ry * ry))


def subplot_title(stem: str) -> str:
    for method in ["godunov", "muscl", "weno2"]:
        if stem.endswith(method) or f"_{method}_" in stem or f"_{method}" in stem:
            return method
    return stem


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", type=Path, nargs="+")
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument("--adv-vel", type=float, default=0.73)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--y-min", type=float, default=0.0)
    parser.add_argument("--y-max", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=360)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(7.5, 5.625))
    parser.add_argument("--total-figsize", type=float, nargs=2, metavar=("W", "H"), default=None)
    parser.add_argument("--interpolation", default="bilinear")
    parser.add_argument("--title", default="2D Gaussian density absolute error")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--title-fontsize", type=int, default=11)
    parser.add_argument("--label-fontsize", type=int, default=10)
    parser.add_argument("--tick-fontsize", type=int, default=8)
    parser.add_argument("--colorbar-fontsize", type=int, default=12)
    args = parser.parse_args()

    for path in args.csv_paths:
        if not path.exists():
            raise SystemExit(f"file not found: {path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    error_grids: list[tuple[Path, list[list[float]]]] = []
    vmax = 0.0
    for path in args.csv_paths:
        rho, nx, ny = load_rho(path)
        err = [[0.0 for _ in range(nx)] for _ in range(ny)]
        for j in range(ny):
            for i in range(nx):
                exact = exact_rho_2d(
                    i,
                    j,
                    nx,
                    ny,
                    args.final_time,
                    args.adv_vel,
                    args.x_min,
                    args.x_max,
                    args.y_min,
                    args.y_max,
                )
                err[j][i] = abs(rho[j][i] - exact)
                vmax = max(vmax, err[j][i])
        error_grids.append((path, err))

    if args.total_figsize is not None:
        canvas_size = tuple(args.total_figsize)
    else:
        canvas_size = (args.figsize[0] * len(error_grids), args.figsize[1])

    fig, axes = plt.subplots(
        1,
        len(error_grids),
        figsize=canvas_size,
        dpi=args.dpi,
        constrained_layout=False,
    )
    if len(error_grids) == 1:
        axes = [axes]

    image = None
    for ax, (path, err) in zip(axes, error_grids):
        image = ax.imshow(
            err,
            origin="lower",
            cmap="viridis",
            aspect="equal",
            interpolation=args.interpolation,
            vmin=0.0,
            vmax=vmax,
        )
        ax.set_title(f"{subplot_title(path.stem)} |rho - rho_exact|", fontsize=args.title_fontsize)
        ax.set_xlabel("i", fontsize=args.label_fontsize)
        ax.set_ylabel("j", fontsize=args.label_fontsize)
        ax.tick_params(axis="both", labelsize=args.tick_fontsize)

    fig.suptitle(args.title, fontsize=args.title_fontsize)
    fig.subplots_adjust(left=0.05, right=0.91, bottom=0.10, top=0.86, wspace=0.18)
    fig.canvas.draw()
    positions = [ax.get_position() for ax in axes]
    y0 = min(position.y0 for position in positions)
    y1 = max(position.y1 for position in positions)
    cax = fig.add_axes([0.925, y0, 0.014, y1 - y0])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label("|rho - rho_exact|", rotation=90, labelpad=10, fontsize=args.colorbar_fontsize)
    colorbar.ax.tick_params(labelsize=args.tick_fontsize)

    if args.output is None:
        root = Path(__file__).resolve().parent.parent
        stem = "__vs__".join(path.stem for path in args.csv_paths)
        output_path = root / "results" / "stageD" / f"{stem}.rho_exact_error.png"
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    print(f"saved exact-error plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
