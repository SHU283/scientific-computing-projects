#!/usr/bin/env python3
"""
Plot two 2D CSV field outputs side by side using a shared color scale.

Usage:
python3 scripts/compare_fields_2d.py \
  results/quadrant_riemann_2d_godunov.csv \
  results/quadrant_riemann_2d_muscl.csv \
  --field rho
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


def min_max(grid: list[list[float]]) -> tuple[float, float]:
    values = [value for row in grid for value in row]
    return min(values), max(values)


def default_stage_d_output(csv_a: Path, csv_b: Path, field: str) -> Path:
    root = Path(__file__).resolve().parent.parent
    return root / "results" / "stageD" / f"{csv_a.stem}__vs__{csv_b.stem}.{field}.png"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_a", type=Path)
    parser.add_argument("csv_b", type=Path)
    parser.add_argument("--field", default="rho", choices=["rho", "u", "v", "p", "energy"])
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(12.0, 5.0))
    parser.add_argument("--interpolation", default="nearest")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    if not args.csv_a.exists():
        raise SystemExit(f"file not found: {args.csv_a}")
    if not args.csv_b.exists():
        raise SystemExit(f"file not found: {args.csv_b}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Install it with: pip install matplotlib")

    grid_a, _, _ = load_field(args.csv_a, args.field)
    grid_b, _, _ = load_field(args.csv_b, args.field)

    min_a, max_a = min_max(grid_a)
    min_b, max_b = min_max(grid_b)
    vmin = min(min_a, min_b)
    vmax = max(max_a, max_b)

    fig, axes = plt.subplots(1, 2, figsize=tuple(args.figsize), dpi=args.dpi, constrained_layout=False)

    image_a = axes[0].imshow(
        grid_a,
        origin="lower",
        cmap="viridis",
        aspect="equal",
        interpolation=args.interpolation,
        vmin=vmin,
        vmax=vmax,
    )
    axes[0].set_title(f"{args.csv_a.name} [{args.field}]")
    axes[0].set_xlabel("i")
    axes[0].set_ylabel("j")

    image_b = axes[1].imshow(
        grid_b,
        origin="lower",
        cmap="viridis",
        aspect="equal",
        interpolation=args.interpolation,
        vmin=vmin,
        vmax=vmax,
    )
    axes[1].set_title(f"{args.csv_b.name} [{args.field}]")
    axes[1].set_xlabel("i")
    axes[1].set_ylabel("j")

    fig.subplots_adjust(left=0.07, right=0.88, bottom=0.11, top=0.90, wspace=0.22)
    cax = fig.add_axes([0.90, 0.16, 0.018, 0.68])
    colorbar = fig.colorbar(image_b, cax=cax)
    colorbar.set_label(args.field, rotation=90, labelpad=10)

    if args.output is None:
        output_path = default_stage_d_output(args.csv_a, args.csv_b, args.field)
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved comparison plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
