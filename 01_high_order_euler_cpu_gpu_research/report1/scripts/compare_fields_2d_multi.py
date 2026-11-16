#!/usr/bin/env python3
"""
Plot multiple 2D CSV fields side by side with a shared colorbar.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

BACKEND_TITLES = {"serial", "openmp", "cuda"}
METHOD_TITLES = {"godunov", "muscl", "weno2", "weno3", "ader3"}


def load_field(path: Path, field_name: str) -> tuple[list[list[float]], int, int]:
    rows: list[tuple[int, int, float]] = []
    max_i = -1
    max_j = -1
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            i = int(row["i"])
            j = int(row["j"])
            value = float(row[field_name])
            rows.append((i, j, value))
            max_i = max(max_i, i)
            max_j = max(max_j, j)
    nx = max_i + 1
    ny = max_j + 1
    grid = [[0.0 for _ in range(nx)] for _ in range(ny)]
    for i, j, value in rows:
        grid[j][i] = value
    return grid, nx, ny


def min_max(grid: list[list[float]]) -> tuple[float, float]:
    values = [value for row in grid for value in row]
    return min(values), max(values)


def subplot_title(stem: str) -> str:
    tokens = stem.split("_")
    for token in tokens:
        if token in BACKEND_TITLES:
            return token
    for token in reversed(tokens):
        if token in METHOD_TITLES:
            return token
    return stem


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", nargs="+", type=Path)
    parser.add_argument("--field", default="rho", choices=["rho", "u", "v", "p", "energy"])
    parser.add_argument("--suptitle", default=None)
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(7.5, 5.625))
    parser.add_argument("--total-figsize", type=float, nargs=2, metavar=("W", "H"), default=None)
    parser.add_argument("--interpolation", default="bilinear")
    parser.add_argument("--title-fontsize", type=int, default=14)
    parser.add_argument("--label-fontsize", type=int, default=12)
    parser.add_argument("--tick-fontsize", type=int, default=10)
    parser.add_argument("--colorbar-fontsize", type=int, default=12)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    if len(args.csv_paths) < 2:
        raise SystemExit("need at least two csv files")
    for path in args.csv_paths:
        if not path.exists():
            raise SystemExit(f"file not found: {path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    grids = []
    vmin = None
    vmax = None
    for path in args.csv_paths:
        grid, _, _ = load_field(path, args.field)
        grids.append((path, grid))
        local_min, local_max = min_max(grid)
        vmin = local_min if vmin is None else min(vmin, local_min)
        vmax = local_max if vmax is None else max(vmax, local_max)

    if args.total_figsize is not None:
        canvas_size = tuple(args.total_figsize)
    else:
        canvas_size = (args.figsize[0] * len(grids), args.figsize[1])

    fig, axes = plt.subplots(
        1,
        len(grids),
        figsize=canvas_size,
        dpi=args.dpi,
        constrained_layout=False,
    )
    if len(grids) == 1:
        axes = [axes]

    image = None
    for ax, (path, grid) in zip(axes, grids):
        image = ax.imshow(
            grid,
            origin="lower",
            cmap="viridis",
            aspect="equal",
            interpolation=args.interpolation,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"{subplot_title(path.stem)} [{args.field}]", fontsize=args.title_fontsize)
        ax.set_xlabel("i", fontsize=args.label_fontsize)
        ax.set_ylabel("j", fontsize=args.label_fontsize)
        ax.tick_params(axis="both", labelsize=args.tick_fontsize)

    top_margin = 0.90
    if args.suptitle:
        fig.suptitle(args.suptitle, fontsize=max(args.title_fontsize + 2, 16))
        top_margin = 0.84

    fig.subplots_adjust(left=0.05, right=0.90, bottom=0.10, top=top_margin, wspace=0.20)
    fig.canvas.draw()
    positions = [ax.get_position() for ax in axes]
    y0 = min(position.y0 for position in positions)
    y1 = max(position.y1 for position in positions)
    cax = fig.add_axes([0.92, y0, 0.015, y1 - y0])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label(args.field, rotation=90, labelpad=10, fontsize=args.colorbar_fontsize)
    colorbar.ax.tick_params(labelsize=args.tick_fontsize)

    if args.output is None:
        root = Path(__file__).resolve().parent.parent
        stem = "__vs__".join(path.stem for path in args.csv_paths)
        output_path = root / "results" / "stageD" / f"{stem}.{args.field}.png"
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    print(f"saved comparison plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

