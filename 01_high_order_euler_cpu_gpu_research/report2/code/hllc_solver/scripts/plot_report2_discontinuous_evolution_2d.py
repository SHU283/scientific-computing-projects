#!/usr/bin/env python3
"""Plot initial, intermediate, and final HLLC discontinuous fields."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CASES = {
    "quadrant_riemann_2d": "Quadrant Riemann density evolution",
    "shock_bubble_2d": "Shock-bubble density evolution",
}


def load_density(path: Path) -> tuple[list[list[float]], int, int]:
    values: list[tuple[int, int, float]] = []
    max_i = -1
    max_j = -1
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            i = int(row["i"])
            j = int(row["j"])
            values.append((i, j, float(row["rho"])))
            max_i = max(max_i, i)
            max_j = max(max_j, j)

    nx = max_i + 1
    ny = max_j + 1
    field = [[0.0 for _ in range(nx)] for _ in range(ny)]
    for i, j, value in values:
        field[j][i] = value
    return field, nx, ny


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=600)
    parser.add_argument(
        "--snapshot-dir",
        type=Path,
        default=Path("results/report2_hllc/discontinuous_evolution_2d/fields"),
    )
    parser.add_argument(
        "--final-dir",
        type=Path,
        default=Path("results/report2_hllc/high_resolution_fields/2d"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/report2_hllc/discontinuous_evolution_2d/figures"),
    )
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = "Arial"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for problem, title in CASES.items():
        paths = [
            args.snapshot_dir
            / f"{problem}_weno3_cuda_n{args.resolution}_t000.csv",
            args.snapshot_dir
            / f"{problem}_weno3_cuda_n{args.resolution}_t010.csv",
            args.final_dir / f"{problem}_weno3_n{args.resolution}.csv",
        ]
        for path in paths:
            if not path.is_file():
                raise SystemExit(f"missing evolution field: {path}")

        fields = [load_density(path)[0] for path in paths]
        vmin = min(value for field in fields for row in field for value in row)
        vmax = max(value for field in fields for row in field for value in row)

        fig, axes = plt.subplots(1, 3, figsize=(16.8, 6.0), dpi=300)
        image = None
        for ax, field, time in zip(axes, fields, (0.0, 0.1, 0.2)):
            image = ax.imshow(
                field,
                origin="lower",
                extent=(0.0, 1.0, 0.0, 1.0),
                cmap="viridis",
                aspect="equal",
                interpolation="nearest",
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_title(f"t = {time:.1f}", fontsize=22, pad=10)
            ax.set_xlabel("x", fontsize=18, labelpad=7)
            ax.set_ylabel("y", fontsize=18, labelpad=7)
            ax.tick_params(axis="both", labelsize=15)

        fig.suptitle(
            f"{title}, WENO-3, {args.resolution} x {args.resolution}",
            fontsize=24,
            y=0.96,
        )
        fig.subplots_adjust(left=0.055, right=0.90, bottom=0.14, top=0.82, wspace=0.22)
        colorbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.025)
        colorbar.set_label("density", fontsize=18, labelpad=10)
        colorbar.ax.tick_params(labelsize=15)

        output = (
            args.output_dir
            / f"{problem}_weno3_evolution_n{args.resolution}.png"
        )
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)
        print(f"saved: {output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
