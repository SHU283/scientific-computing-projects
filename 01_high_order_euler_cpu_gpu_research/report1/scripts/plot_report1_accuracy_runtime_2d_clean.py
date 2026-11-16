#!/usr/bin/env python3
"""
Cleaner Report 1 accuracy-runtime plot for the 2D Gaussian benchmark.

The plot is split by numerical method so the extra n=25 point can be shown
without crowding every backend/method label into a single panel.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METHODS = ["godunov", "muscl", "weno2"]
BACKENDS = ["serial", "openmp", "cuda"]


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def pretty_metric(metric: str) -> str:
    return {
        "l1_rho": "L1 density error",
        "l2_rho": "L2 density error",
        "linf_rho": "Linf density error",
        "best_seconds": "best runtime (s)",
        "mean_seconds": "mean runtime (s)",
        "median_seconds": "median runtime (s)",
    }.get(metric, metric)


def make_method_axes(plt, methods: list[str], figsize: tuple[float, float], dpi: int):
    """Use a 2+1 method layout so each panel stays readable in the report."""
    fig = plt.figure(figsize=figsize, dpi=dpi)
    grid = fig.add_gridspec(2, 4)
    slots = [
        grid[0, 0:2],
        grid[0, 2:4],
        grid[1, 1:3],
    ]
    axes = []
    for index in range(len(methods)):
        if index == 0:
            axes.append(fig.add_subplot(slots[index]))
        else:
            axes.append(fig.add_subplot(slots[index], sharey=axes[0]))
    return fig, axes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--error-metric", default="l1_rho", choices=["l1_rho", "l2_rho", "linf_rho"])
    parser.add_argument("--time-metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(11.0, 8.0))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--no-n-labels", action="store_true")
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    rows = [row for row in load_rows(args.csv_path) if row["problem"] == "gaussian_advection_2d"]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["method"], row["backend"])].append(row)

    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    markers = {"serial": "o", "openmp": "s", "cuda": "^"}

    fig, axes = make_method_axes(plt, METHODS, tuple(args.figsize), args.dpi)
    for ax, method in zip(axes, METHODS):
        for backend in BACKENDS:
            series = grouped.get((method, backend), [])
            if not series:
                continue
            series = sorted(series, key=lambda row: int(row["nx"]))
            xs = [float(row[args.time_metric]) for row in series]
            ys = [float(row[args.error_metric]) for row in series]
            ns = [int(row["nx"]) for row in series]
            ax.plot(
                xs,
                ys,
                marker=markers[backend],
                markersize=6,
                linewidth=1.8,
                color=colors[backend],
                label=backend,
            )
            if not args.no_n_labels:
                for x, y, n in zip(xs, ys, ns):
                    ax.annotate(f"n={n}", (x, y), textcoords="offset points", xytext=(5, 4), fontsize=8)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(method, fontsize=14)
        ax.set_xlabel(pretty_metric(args.time_metric), fontsize=12)
        ax.grid(True, alpha=0.25)
        ax.tick_params(axis="both", labelsize=10)

    axes[0].set_ylabel(pretty_metric(args.error_metric), fontsize=12)
    axes[-1].legend(title="backend", fontsize=10, title_fontsize=11, loc="best")
    fig.suptitle("2D Gaussian accuracy-runtime tradeoff", fontsize=16)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.91, wspace=0.40, hspace=0.48)

    if args.output is None:
        output_path = args.csv_path.with_suffix(f".clean.{args.time_metric}.{args.error_metric}.png")
    else:
        output_path = args.output

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved clean accuracy-runtime plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
