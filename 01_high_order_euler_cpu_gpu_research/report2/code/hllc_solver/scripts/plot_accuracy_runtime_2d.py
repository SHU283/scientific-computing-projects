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


def pretty_metric(metric: str) -> str:
    return {
        "l1_rho": "L1 density error",
        "l2_rho": "L2 density error",
        "linf_rho": "Linf density error",
        "best_seconds": "best runtime (s)",
        "mean_seconds": "mean runtime (s)",
        "median_seconds": "median runtime (s)",
    }.get(metric, metric)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--error-metric", default="l1_rho", choices=["l1_rho", "l2_rho", "linf_rho"])
    parser.add_argument("--time-metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 6.0))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--no-point-labels", action="store_true")
    parser.add_argument("--label-mode", default="all", choices=["all", "serial-only"])
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
    markers = {"godunov": "o", "muscl": "s", "weno2": "^"}
    label_offsets = {
        ("serial", "godunov", 50): (4, -8),
        ("openmp", "godunov", 50): (-20, 8),
        ("cuda", "godunov", 50): (4, 8),
        ("serial", "godunov", 100): (4, 6),
        ("openmp", "godunov", 100): (-24, 6),
        ("cuda", "godunov", 100): (4, -9),
        ("serial", "godunov", 200): (4, 5),
        ("openmp", "godunov", 200): (-10, 14),
        ("cuda", "godunov", 200): (8, -14),
        ("serial", "muscl", 50): (4, 8),
        ("openmp", "muscl", 50): (-22, 8),
        ("cuda", "muscl", 50): (8, 14),
        ("serial", "muscl", 100): (-31, -10),
        ("openmp", "muscl", 100): (6, -12),
        ("cuda", "muscl", 100): (5, -18),
        ("serial", "muscl", 200): (4, 6),
        ("openmp", "muscl", 200): (4, -10),
        ("cuda", "muscl", 200): (-28, 9),
        ("serial", "weno2", 50): (4, -10),
        ("openmp", "weno2", 50): (-22, -12),
        ("cuda", "weno2", 50): (4, 6),
        ("serial", "weno2", 100): (4, 10),
        ("openmp", "weno2", 100): (-28, 8),
        ("cuda", "weno2", 100): (4, -11),
        ("serial", "weno2", 200): (4, 5),
        ("openmp", "weno2", 200): (4, -12),
        ("cuda", "weno2", 200): (-26, 9),
    }

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)

    for row in rows:
        x = float(row[args.time_metric])
        y = float(row[args.error_metric])
        backend = row["backend"]
        method = row["method"]
        n = int(row["nx"])

        ax.scatter(
            x,
            y,
            s=55,
            color=colors.get(backend, "#333333"),
            marker=markers.get(method, "o"),
            alpha=0.9,
        )
        show_label = not args.no_point_labels
        if args.label_mode == "serial-only" and backend != "serial":
            show_label = False

        if show_label:
            offset = label_offsets.get((backend, method, n), (6, 4))
            ax.annotate(
                f"n={n}",
                (x, y),
                textcoords="offset points",
                xytext=offset,
                fontsize=6.5,
                color="#333333",
                clip_on=False,
            )

    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.set_xlabel(pretty_metric(args.time_metric))
    ax.set_ylabel(pretty_metric(args.error_metric))
    ax.set_title("Accuracy vs runtime: 2D Gaussian advection")
    ax.grid(True, alpha=0.25)
    ax.margins(x=0.18, y=0.18)

    from matplotlib.lines import Line2D

    backend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[name], markersize=8, label=name)
        for name in ["serial", "openmp", "cuda"]
    ]
    method_handles = [
        Line2D([0], [0], marker=markers[name], color="#444444", linestyle="None", markersize=8, label=name)
        for name in ["godunov", "muscl", "weno2"]
    ]
    legend1 = ax.legend(handles=backend_handles, title="backend", loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.add_artist(legend1)
    ax.legend(handles=method_handles, title="method", loc="upper left", bbox_to_anchor=(1.02, 0.62))

    fig.tight_layout(rect=(0, 0, 0.82, 1))

    if args.output is None:
        output_path = args.csv_path.with_suffix(f".{args.time_metric}.{args.error_metric}.png")
    else:
        output_path = args.output

    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
