#!/usr/bin/env python3
"""
Read a 2D CSV field and render a density plot.

Usage:
python3 scripts/plot_field_2d.py results/shock_bubble_2d_godunov.csv
python3 scripts/plot_field_2d.py results/quadrant_riemann_2d_muscl.csv --field=rho
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--field", default="rho", choices=["rho", "u", "v", "p", "energy"])
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=720)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(7.5, 5.625))
    parser.add_argument("--interpolation", default="bilinear")
    parser.add_argument("--title-fontsize", type=int, default=16)
    parser.add_argument("--label-fontsize", type=int, default=12)
    parser.add_argument("--tick-fontsize", type=int, default=10)
    parser.add_argument("--colorbar-fontsize", type=int, default=12)
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required; you can first run: pip install matplotlib")

    grid, nx, ny = load_field(args.csv_path, args.field)

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)
    image = ax.imshow(
        grid,
        origin="lower",
        cmap="viridis",
        aspect="equal",
        interpolation=args.interpolation,
    )
    ax.set_title(f"{args.csv_path.name} [{args.field}]", fontsize=args.title_fontsize)
    ax.set_xlabel("i", fontsize=args.label_fontsize)
    ax.set_ylabel("j", fontsize=args.label_fontsize)
    ax.tick_params(axis="both", labelsize=args.tick_fontsize)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.9)
    colorbar.set_label(args.field, fontsize=args.colorbar_fontsize)
    colorbar.ax.tick_params(labelsize=args.tick_fontsize)
    fig.tight_layout()

    if args.output is None:
        output_path = args.csv_path.with_suffix(f".{args.field}.png")
    else:
        output_path = args.output

    fig.savefig(output_path)
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

