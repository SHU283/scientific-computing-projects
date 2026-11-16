#!/usr/bin/env python3
"""
Read a 2D discontinuous runtime CSV and plot:
1. runtime comparison
2. speedup comparison
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def pretty_label(text: str) -> str:
    return text.replace("_", " ")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(13.5, 4.8))
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.ticker import FuncFormatter
        from matplotlib.patches import Patch
    except ImportError:
        raise SystemExit("matplotlib and numpy are required. Try: pip install matplotlib numpy")

    plt.rcParams["text.parse_math"] = False

    rows = load_rows(args.csv_path)
    if not rows:
        raise SystemExit("runtime csv is empty")

    cases = []
    methods = []
    backends = ["serial", "openmp", "cuda"]
    for row in rows:
        if row["case"] not in cases:
            cases.append(row["case"])
        if row["method"] not in methods:
            methods.append(row["method"])

    values: dict[tuple[str, str, str], float] = {}
    serial_baseline: dict[tuple[str, str], float] = {}
    for row in rows:
        key = (row["case"], row["method"], row["backend"])
        value = float(row[args.metric])
        values[key] = value
        if row["backend"] == "serial":
            serial_baseline[(row["case"], row["method"])] = value

    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    case_labels = {
        "quadrant_riemann_2d": "quadrant riemann 2d",
        "shock_bubble_2d": "shock bubble 2d",
    }
    width = 0.22

    runtime_fig, runtime_axes = plt.subplots(1, len(methods), figsize=tuple(args.figsize), dpi=args.dpi, sharey=True)
    if len(methods) == 1:
        runtime_axes = [runtime_axes]

    x = np.arange(len(cases))
    for ax, method in zip(runtime_axes, methods):
        for idx, backend in enumerate(backends):
            ys = [values.get((case, method, backend), 0.0) for case in cases]
            ax.bar(x + (idx - 1) * width, ys, width=width, color=colors[backend], label=backend)

        ax.set_title(f"Runtime: {method}")
        ax.set_xticks(x)
        ax.set_xticklabels([case_labels.get(case, pretty_label(case)) for case in cases], rotation=0)
        ax.set_yscale("log", base=10)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
        ax.grid(True, axis="y", alpha=0.25)

    runtime_axes[0].set_ylabel(f"{pretty_label(args.metric)} (s)")
    runtime_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    runtime_fig.legend(handles=runtime_handles, loc="upper center", ncol=len(backends), frameon=True)

    runtime_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.82, wspace=0.18)
    runtime_out = args.csv_path.with_suffix(f".{args.metric}.runtime.png")
    runtime_fig.savefig(runtime_out)
    print(f"saved plot: {runtime_out}")

    speedup_fig, speedup_axes = plt.subplots(1, len(methods), figsize=tuple(args.figsize), dpi=args.dpi, sharey=True)
    if len(methods) == 1:
        speedup_axes = [speedup_axes]

    for ax, method in zip(speedup_axes, methods):
        for idx, backend in enumerate(backends):
            ys: list[float] = []
            for case in cases:
                baseline = serial_baseline.get((case, method))
                value = values.get((case, method, backend))
                if baseline is None or value is None or value <= 0.0:
                    ys.append(0.0)
                else:
                    ys.append(baseline / value)
            ax.bar(x + (idx - 1) * width, ys, width=width, color=colors[backend], label=backend)

        ax.set_title(f"Speedup: {method}")
        ax.set_xticks(x)
        ax.set_xticklabels([case_labels.get(case, pretty_label(case)) for case in cases], rotation=0)
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
        ax.grid(True, axis="y", alpha=0.25)

    speedup_axes[0].set_ylabel("speedup vs serial")
    speedup_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    speedup_fig.legend(handles=speedup_handles, loc="upper center", ncol=len(backends), frameon=True)

    speedup_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.82, wspace=0.18)
    speedup_out = args.csv_path.with_suffix(f".{args.metric}.speedup.png")
    speedup_fig.savefig(speedup_out)
    print(f"saved plot: {speedup_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

