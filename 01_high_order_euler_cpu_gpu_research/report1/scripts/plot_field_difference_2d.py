#!/usr/bin/env python3
"""
Plot the difference between two 2D field CSV files.

Default behavior uses absolute difference with the same viridis colormap
as the side-by-side comparison plots, so the visual style stays consistent.
Use --signed to switch back to a signed difference plot with a diverging map.
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


def default_stage_d_output(csv_a: Path, csv_b: Path, field: str) -> Path:
    root = Path(__file__).resolve().parent.parent
    return root / "results" / "stageD" / f"{csv_a.stem}__to__{csv_b.stem}.diff_{field}.png"


def choose_title_fontsize(title: str, requested: int) -> int:
    title_len = len(title)
    if title_len <= 60:
        return requested
    if title_len <= 78:
        return max(12, requested - 1)
    if title_len <= 96:
        return max(11, requested - 2)
    return max(10, requested - 3)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_a", type=Path)
    parser.add_argument("csv_b", type=Path)
    parser.add_argument("--field", default="rho", choices=["rho", "u", "v", "p", "energy"])
    parser.add_argument("--signed", action="store_true", help="plot signed difference with a diverging colormap")
    parser.add_argument("--cmap", default=None, help="override colormap name")
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(12.0, 9.0))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--interpolation", default="bilinear")
    parser.add_argument("--title-fontsize", type=int, default=14)
    parser.add_argument("--label-fontsize", type=int, default=12)
    parser.add_argument("--tick-fontsize", type=int, default=10)
    parser.add_argument("--colorbar-fontsize", type=int, default=12)
    args = parser.parse_args()

    if not args.csv_a.exists():
        raise SystemExit(f"file not found: {args.csv_a}")
    if not args.csv_b.exists():
        raise SystemExit(f"file not found: {args.csv_b}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    grid_a, nx_a, ny_a = load_field(args.csv_a, args.field)
    grid_b, nx_b, ny_b = load_field(args.csv_b, args.field)
    if nx_a != nx_b or ny_a != ny_b:
        raise SystemExit("grid size mismatch between the two csv files")

    raw_diff = [[grid_b[j][i] - grid_a[j][i] for i in range(nx_a)] for j in range(ny_a)]
    max_abs = max(abs(value) for row in raw_diff for value in row)

    if args.signed:
        display_grid = raw_diff
        vmin = -(max_abs if max_abs > 0.0 else 1.0)
        vmax = max_abs if max_abs > 0.0 else 1.0
        cmap = args.cmap or "coolwarm"
        title = f"{args.csv_b.name} - {args.csv_a.name} [{args.field}]"
        colorbar_label = f"delta {args.field}"
    else:
        display_grid = [[abs(value) for value in row] for row in raw_diff]
        vmin = 0.0
        vmax = max_abs if max_abs > 0.0 else 1.0
        cmap = args.cmap or "viridis"
        title = f"|{args.csv_b.name} - {args.csv_a.name}| [{args.field}]"
        colorbar_label = f"|delta {args.field}|"

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi, constrained_layout=False)
    image = ax.imshow(
        display_grid,
        origin="lower",
        cmap=cmap,
        aspect="equal",
        interpolation=args.interpolation,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title, fontsize=choose_title_fontsize(title, args.title_fontsize))
    ax.set_xlabel("i", fontsize=args.label_fontsize)
    ax.set_ylabel("j", fontsize=args.label_fontsize)
    ax.tick_params(axis="both", labelsize=args.tick_fontsize)

    fig.subplots_adjust(left=0.10, right=0.84, bottom=0.11, top=0.90)
    fig.canvas.draw()
    position = ax.get_position()
    cax = fig.add_axes([0.85, position.y0, 0.023, position.y1 - position.y0])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label(colorbar_label, rotation=90, labelpad=10, fontsize=args.colorbar_fontsize)
    colorbar.ax.tick_params(labelsize=args.tick_fontsize)

    if args.output is None:
        output_path = default_stage_d_output(args.csv_a, args.csv_b, args.field)
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    print(f"saved difference plot: {output_path}")
    print(f"max |delta {args.field}| = {max_abs:.12e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

