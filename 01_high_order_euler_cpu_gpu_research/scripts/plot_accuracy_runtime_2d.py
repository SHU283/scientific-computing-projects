#!/usr/bin/env python3
"""
Read accuracy_runtime_2d.csv and plot error vs runtime.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--error-metric", default="l1_rho", choices=["l1_rho", "l2_rho", "linf_rho"])
    parser.add_argument("--time-metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 6.0))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    plt.rcParams["mathtext.fontset"] = "dejavusans"
    plt.rcParams["font.family"] = "DejaVu Sans"

    rows = load_rows(args.csv_path)
    if not rows:
        raise SystemExit("accuracy-runtime csv is empty")

    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    markers = {"godunov": "o", "muscl": "s"}

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)

    for row in rows:
        x = float(row[args.time_metric])
        y = float(row[args.error_metric])
        backend = row["backend"]
        method = row["method"]
        n = row["nx"]

        ax.scatter(
            x,
            y,
            s=80,
            color=colors.get(backend, "#333333"),
            marker=markers.get(method, "o"),
            alpha=0.9,
        )
        ax.annotate(
            f"{backend}-{method}-n{n}",
            (x, y),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=8,
        )

    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.set_xlabel(f"{args.time_metric} (s)")
    ax.set_ylabel(args.error_metric)
    ax.set_title("2D Accuracy vs Runtime")
    ax.grid(True, alpha=0.25)

    from matplotlib.lines import Line2D

    backend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[name], markersize=8, label=name)
        for name in ["serial", "openmp", "cuda"]
    ]
    method_handles = [
        Line2D([0], [0], marker=markers[name], color="#444444", linestyle="None", markersize=8, label=name)
        for name in ["godunov", "muscl"]
    ]
    legend1 = ax.legend(handles=backend_handles, title="backend", loc="upper right")
    ax.add_artist(legend1)
    ax.legend(handles=method_handles, title="method", loc="lower left")

    fig.tight_layout()

    if args.output is None:
        output_path = args.csv_path.with_suffix(f".{args.time_metric}.{args.error_metric}.png")
    else:
        output_path = args.output

    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
