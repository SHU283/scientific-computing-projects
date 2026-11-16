#!/usr/bin/env python3
"""
Read a 2D runtime scaling CSV and plot:
1. runtime vs resolution
2. speedup vs resolution
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def format_runtime_tick(value: float) -> str:
    if value <= 0.0:
        return ""
    if value >= 1.0:
        return f"{value:g}"
    if value >= 1e-2:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    if value >= 1e-3:
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return f"{value:.1e}"


def pretty_metric(metric: str) -> str:
    return {
        "best_seconds": "best runtime (s)",
        "mean_seconds": "mean runtime (s)",
        "median_seconds": "median runtime (s)",
    }.get(metric, metric)


def make_method_axes(plt, methods: list[str], figsize: tuple[float, float], dpi: int):
    """Use a larger 2+1 layout so each method panel remains readable in print."""
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
    parser.add_argument("--metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(11.0, 8.0))
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        from matplotlib.ticker import FuncFormatter, NullFormatter
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    plt.rcParams["text.parse_math"] = False

    rows = load_rows(args.csv_path)
    if not rows:
        raise SystemExit("runtime scaling csv is empty")

    by_method_backend: dict[tuple[str, str], list[tuple[int, float]]] = defaultdict(list)
    serial_baseline: dict[tuple[str, int], float] = {}

    for row in rows:
        method = row["method"]
        backend = row["backend"]
        nx = int(row["nx"])
        value = float(row[args.metric])
        by_method_backend[(method, backend)].append((nx, value))
        if backend == "serial":
            serial_baseline[(method, nx)] = value

    for key in by_method_backend:
        by_method_backend[key].sort(key=lambda item: item[0])

    methods = sorted({row["method"] for row in rows})
    resolutions = sorted({int(row["nx"]) for row in rows})
    backends = ["serial", "openmp", "cuda"]
    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    markers = {"serial": "o", "openmp": "s", "cuda": "^"}

    runtime_fig, runtime_axes = make_method_axes(plt, methods, tuple(args.figsize), args.dpi)

    for ax, method in zip(runtime_axes, methods):
        for backend in backends:
            points = by_method_backend.get((method, backend), [])
            if not points:
                continue
            xs = [nx for nx, _ in points]
            ys = [val for _, val in points]
            ax.plot(xs, ys, marker=markers[backend], color=colors[backend], linewidth=2, label=backend)

        ax.set_title(f"Runtime Scaling: {method}")
        ax.set_xlabel("resolution n")
        ax.set_ylabel(pretty_metric(args.metric))
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=10)
        ax.set_xticks(resolutions)
        ax.set_yticks([0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):d}" if int(round(x)) in resolutions else ""))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: format_runtime_tick(y)))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", which="both", labelleft=True)
        ax.grid(True, alpha=0.25)
    runtime_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    runtime_fig.legend(handles=runtime_handles, loc="upper center", ncol=len(backends), frameon=True)

    runtime_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.91, wspace=0.40, hspace=0.48)
    runtime_out = args.csv_path.with_suffix(f".{args.metric}.runtime.png")
    runtime_fig.savefig(runtime_out, bbox_inches="tight")
    print(f"saved plot: {runtime_out}")

    speedup_fig, speedup_axes = make_method_axes(plt, methods, tuple(args.figsize), args.dpi)

    for ax, method in zip(speedup_axes, methods):
        for backend in backends:
            points = by_method_backend.get((method, backend), [])
            if not points:
                continue
            xs: list[int] = []
            ys: list[float] = []
            for nx, val in points:
                baseline = serial_baseline.get((method, nx))
                if baseline is None or val <= 0.0:
                    continue
                xs.append(nx)
                ys.append(baseline / val)

            ax.plot(xs, ys, marker=markers[backend], color=colors[backend], linewidth=2, label=backend)

        ax.set_title(f"Speedup Scaling: {method}")
        ax.set_xlabel("resolution n")
        ax.set_ylabel("speedup vs serial")
        ax.set_xscale("log", base=2)
        ax.set_xticks(resolutions)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):d}" if int(round(x)) in resolutions else ""))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", which="both", labelleft=True)
        ax.grid(True, alpha=0.25)
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
    speedup_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    speedup_fig.legend(handles=speedup_handles, loc="upper center", ncol=len(backends), frameon=True)

    speedup_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.91, wspace=0.40, hspace=0.48)
    speedup_out = args.csv_path.with_suffix(f".{args.metric}.speedup.png")
    speedup_fig.savefig(speedup_out, bbox_inches="tight")
    print(f"saved plot: {speedup_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
