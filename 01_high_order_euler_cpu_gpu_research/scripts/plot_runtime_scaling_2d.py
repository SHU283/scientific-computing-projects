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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 5.8))
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
    backends = ["serial", "openmp", "cuda"]
    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    markers = {"serial": "o", "openmp": "s", "cuda": "^"}

    runtime_fig, runtime_axes = plt.subplots(
        1,
        len(methods),
        figsize=(args.figsize[0] * len(methods), args.figsize[1]),
        dpi=args.dpi,
        sharey=True,
    )
    if len(methods) == 1:
        runtime_axes = [runtime_axes]

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
        ax.set_ylabel(f"{args.metric} (s)")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=10)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):d}" if x >= 1 else ""))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: format_runtime_tick(y)))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", which="both", labelleft=True)
        ax.grid(True, alpha=0.25)
    runtime_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    runtime_fig.legend(handles=runtime_handles, loc="upper center", ncol=len(backends), frameon=True)

    runtime_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.85, wspace=0.16)
    runtime_out = args.csv_path.with_suffix(f".{args.metric}.runtime.png")
    runtime_fig.savefig(runtime_out, bbox_inches="tight")
    print(f"saved plot: {runtime_out}")

    speedup_fig, speedup_axes = plt.subplots(
        1,
        len(methods),
        figsize=(args.figsize[0] * len(methods), args.figsize[1]),
        dpi=args.dpi,
        sharey=True,
    )
    if len(methods) == 1:
        speedup_axes = [speedup_axes]

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
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):d}" if x >= 1 else ""))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", which="both", labelleft=True)
        ax.grid(True, alpha=0.25)
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
    speedup_handles = [Patch(facecolor=colors[name], label=name) for name in backends]
    speedup_fig.legend(handles=speedup_handles, loc="upper center", ncol=len(backends), frameon=True)

    speedup_fig.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.85, wspace=0.16)
    speedup_out = args.csv_path.with_suffix(f".{args.metric}.speedup.png")
    speedup_fig.savefig(speedup_out, bbox_inches="tight")
    print(f"saved plot: {speedup_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
