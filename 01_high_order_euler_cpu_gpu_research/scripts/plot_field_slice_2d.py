#!/usr/bin/env python3
"""
Extract and plot a horizontal or vertical slice from one or more 2D CSV fields.

Usage:
python3 scripts/plot_field_slice_2d.py \
  results/shock_bubble_2d_godunov.csv \
  results/shock_bubble_2d_muscl.csv \
  --field rho --axis y --index 50
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


def default_stage_d_output(csv_paths: list[Path], axis: str, index: int, field: str) -> Path:
    root = Path(__file__).resolve().parent.parent
    stem = "__vs__".join(path.stem for path in csv_paths)
    return root / "results" / "stageD" / f"{stem}.slice_{axis}{index}_{field}.png"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", type=Path, nargs="+")
    parser.add_argument("--field", default="rho", choices=["rho", "u", "v", "p", "energy"])
    parser.add_argument("--axis", default="y", choices=["x", "y"])
    parser.add_argument("--index", type=int, default=50)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 5.2))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    for path in args.csv_paths:
        if not path.exists():
            raise SystemExit(f"file not found: {path}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Install it with: pip install matplotlib")

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)
    colors = ["#1f4e79", "#d97706", "#0f766e", "#b91c1c", "#6d28d9"]

    for idx_path, path in enumerate(args.csv_paths):
        grid, nx, ny = load_field(path, args.field)

        if args.axis == "y":
            if args.index < 0 or args.index >= ny:
                raise SystemExit(f"slice index out of range for y-axis: {args.index}")
            xs = list(range(nx))
            ys = grid[args.index]
            xlabel = "i"
            title_axis = f"y-slice j={args.index}"
        else:
            if args.index < 0 or args.index >= nx:
                raise SystemExit(f"slice index out of range for x-axis: {args.index}")
            xs = list(range(ny))
            ys = [grid[j][args.index] for j in range(ny)]
            xlabel = "j"
            title_axis = f"x-slice i={args.index}"

        ax.plot(xs, ys, linewidth=2, color=colors[idx_path % len(colors)], label=path.stem)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(args.field)
    ax.set_title(f"{title_axis} [{args.field}]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if args.output is None:
        output_path = default_stage_d_output(args.csv_paths, args.axis, args.index, args.field)
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved slice plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
