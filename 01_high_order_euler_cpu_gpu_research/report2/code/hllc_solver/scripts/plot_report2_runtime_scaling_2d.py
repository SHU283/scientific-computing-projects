#!/usr/bin/env python3
"""Plot matched-backend Report 2 runtime and speedup scaling."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
METHOD_LABELS = {
    "weno2": "WENO-2",
    "weno3": "WENO-3",
    "ader3": "ADER-3",
}
BACKENDS = ("serial", "openmp", "cuda")
COLORS = {
    "serial": "#1f4e79",
    "openmp": "#d97706",
    "cuda": "#0f766e",
}
MARKERS = {
    "serial": "o",
    "openmp": "s",
    "cuda": "^",
}
SPEEDUP_LABEL_Y_OFFSETS = {
    ("weno2", "openmp"): (10, -10, 8, -8, 5),
    ("weno2", "cuda"): (6, 8, -8, 8, -8),
    ("weno3", "openmp"): (10, -10, 8, -8, 5),
    ("weno3", "cuda"): (6, 8, -8, 8, -8),
    ("ader3", "openmp"): (-8, 12, -10, 5, 5),
    ("ader3", "cuda"): (6, -14, 8, -8, -8),
}
FIGURE_SIZE = (15.2, 5.6)
FIGURE_TITLE_SIZE = 21
PANEL_TITLE_SIZE = 18
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 15
X_TICK_LABEL_SIZE = 15
LEGEND_SIZE = 15
VALUE_LABEL_SIZE = 11


def log_tick(value: float, _position: int) -> str:
    if value <= 0:
        return ""
    exponent = round(math.log10(value))
    if not math.isclose(value, 10**exponent, rel_tol=1e-8, abs_tol=0.0):
        return ""
    if exponent == 0:
        return "1"
    if exponent == 1:
        return "10"
    if exponent == 2:
        return "100"
    if exponent == -1:
        return "0.1"
    if exponent == -2:
        return "0.01"
    return f"1e{exponent}"


def plain_log_y(axis) -> None:
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    axis.yaxis.set_major_locator(LogLocator(base=10))
    axis.yaxis.set_major_formatter(FuncFormatter(log_tick))
    axis.yaxis.set_minor_formatter(NullFormatter())


def separate_close_x_tick_labels(axis, values: list[int]) -> None:
    """Separate the close 25/50 labels while retaining true grid spacing."""
    if len(values) < 2 or values[1] - values[0] > 25:
        return
    labels = axis.get_xticklabels()
    if len(labels) >= 2:
        labels[0].set_horizontalalignment("right")
        labels[1].set_horizontalalignment("center")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument(
        "--metric",
        default="best_seconds",
        choices=("best_seconds", "mean_seconds", "median_seconds"),
    )
    parser.add_argument("--dpi", type=int, default=240)
    args = parser.parse_args()
    metric_label = args.metric.replace("_seconds", "")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as path_effects
        plt.rcParams["text.parse_math"] = False
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    with args.csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("runtime CSV is empty")

    values: dict[tuple[str, str], list[tuple[int, float]]] = defaultdict(list)
    serial: dict[tuple[str, int], float] = {}
    for row in rows:
        method = row["method"]
        backend = row["backend"]
        resolution = int(row["nx"])
        runtime = float(row[args.metric])
        values[(method, backend)].append((resolution, runtime))
        if backend == "serial":
            serial[(method, resolution)] = runtime
    for points in values.values():
        points.sort()
    all_runtimes = [
        runtime for points in values.values() for _resolution, runtime in points
    ]
    runtime_limits = (min(all_runtimes) * 0.65, max(all_runtimes) * 1.45)

    methods = [
        method
        for method in METHODS
        if any(row["method"] == method for row in rows)
    ]
    resolutions = sorted({int(row["nx"]) for row in rows})
    handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[backend],
            marker=MARKERS[backend],
            linewidth=2,
            label=backend,
        )
        for backend in BACKENDS
    ]

    runtime_fig, runtime_axes = plt.subplots(
        1,
        len(methods),
        figsize=FIGURE_SIZE,
        dpi=args.dpi,
        sharey=True,
    )
    if len(methods) == 1:
        runtime_axes = [runtime_axes]
    for axis, method in zip(runtime_axes, methods):
        for backend in BACKENDS:
            points = values.get((method, backend), [])
            if not points:
                continue
            axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                color=COLORS[backend],
                marker=MARKERS[backend],
                linewidth=2,
                markersize=6,
            )
        axis.set_title(METHOD_LABELS[method], fontsize=PANEL_TITLE_SIZE, pad=10)
        axis.set_yscale("log")
        axis.set_ylim(*runtime_limits)
        plain_log_y(axis)
        axis.set_xticks(resolutions)
        separate_close_x_tick_labels(axis, resolutions)
        axis.set_xlabel(
            "resolution n", fontsize=AXIS_LABEL_SIZE, labelpad=7
        )
        axis.tick_params(axis="x", labelsize=X_TICK_LABEL_SIZE, pad=4)
        axis.tick_params(axis="y", labelsize=TICK_LABEL_SIZE, pad=4)
        for spine in axis.spines.values():
            spine.set_linewidth(1.1)
        axis.grid(True, alpha=0.25)
    runtime_axes[0].set_ylabel(
        f"{metric_label} runtime (s)",
        fontsize=AXIS_LABEL_SIZE,
        labelpad=9,
    )
    runtime_fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=3,
        frameon=True,
        fontsize=LEGEND_SIZE,
    )
    runtime_fig.suptitle(
        "Report 2: 2D Gaussian runtime scaling",
        y=0.985,
        fontsize=FIGURE_TITLE_SIZE,
    )
    runtime_fig.subplots_adjust(
        left=0.075, right=0.985, bottom=0.15, top=0.72, wspace=0.08
    )
    runtime_path = args.csv_path.with_suffix(f".{args.metric}.runtime_report2.png")
    runtime_fig.savefig(runtime_path)
    print(f"saved plot: {runtime_path}")

    speedup_fig, speedup_axes = plt.subplots(
        1,
        len(methods),
        figsize=FIGURE_SIZE,
        dpi=args.dpi,
        sharey=True,
    )
    if len(methods) == 1:
        speedup_axes = [speedup_axes]
    all_speedups: list[float] = []
    for axis, method in zip(speedup_axes, methods):
        for backend in ("openmp", "cuda"):
            points = values.get((method, backend), [])
            xs = []
            ys = []
            for resolution, runtime in points:
                baseline = serial.get((method, resolution))
                if baseline is None:
                    continue
                xs.append(resolution)
                ys.append(baseline / runtime)
            all_speedups.extend(ys)
            axis.plot(
                xs,
                ys,
                color=COLORS[backend],
                marker=MARKERS[backend],
                linewidth=2,
                markersize=6,
                label=backend,
            )
            for point_index, (x, y) in enumerate(zip(xs, ys)):
                offsets = SPEEDUP_LABEL_Y_OFFSETS[(method, backend)]
                y_offset = offsets[min(point_index, len(offsets) - 1)]
                vertical_alignment = "bottom" if y_offset > 0 else "top"

                x_offset = 0
                horizontal_alignment = "center"
                if point_index == 0:
                    x_offset = 4
                    horizontal_alignment = "left"
                elif point_index == len(xs) - 1:
                    x_offset = -4
                    horizontal_alignment = "right"

                annotation = axis.annotate(
                    f"{y:.2f}x" if y < 0.1 else f"{y:.1f}x",
                    (x, y),
                    xytext=(x_offset, y_offset),
                    textcoords="offset points",
                    fontsize=VALUE_LABEL_SIZE,
                    fontweight="normal",
                    color=COLORS[backend],
                    horizontalalignment=horizontal_alignment,
                    verticalalignment=vertical_alignment,
                    clip_on=True,
                )
                annotation.set_path_effects(
                    [path_effects.withStroke(linewidth=1.0, foreground="white")]
                )
        axis.axhline(1.0, color="#666666", linestyle="--", linewidth=1)
        axis.set_title(METHOD_LABELS[method], fontsize=PANEL_TITLE_SIZE, pad=10)
        axis.set_yscale("log")
        plain_log_y(axis)
        axis.set_xticks(resolutions)
        separate_close_x_tick_labels(axis, resolutions)
        axis.set_xlabel(
            "resolution n", fontsize=AXIS_LABEL_SIZE, labelpad=7
        )
        axis.tick_params(axis="x", labelsize=X_TICK_LABEL_SIZE, pad=4)
        axis.tick_params(axis="y", labelsize=TICK_LABEL_SIZE, pad=4)
        for spine in axis.spines.values():
            spine.set_linewidth(1.1)
        axis.grid(True, which="both", alpha=0.25)
    speedup_limits = (min(all_speedups) * 0.65, max(all_speedups) * 1.45)
    for axis in speedup_axes:
        axis.set_ylim(*speedup_limits)
    speedup_axes[0].set_ylabel(
        "speedup versus serial", fontsize=AXIS_LABEL_SIZE, labelpad=9
    )
    speedup_handles = [handles[1], handles[2]]
    speedup_fig.legend(
        handles=speedup_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=2,
        frameon=True,
        fontsize=LEGEND_SIZE,
    )
    speedup_fig.suptitle(
        "Report 2: 2D Gaussian backend speedup",
        y=0.985,
        fontsize=FIGURE_TITLE_SIZE,
    )
    speedup_fig.subplots_adjust(
        left=0.075, right=0.985, bottom=0.15, top=0.72, wspace=0.08
    )
    speedup_path = args.csv_path.with_suffix(
        f".{args.metric}.speedup_report2.png"
    )
    speedup_fig.savefig(speedup_path)
    print(f"saved plot: {speedup_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
