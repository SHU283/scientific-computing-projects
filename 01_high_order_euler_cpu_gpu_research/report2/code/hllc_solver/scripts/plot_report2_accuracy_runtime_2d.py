#!/usr/bin/env python3
"""Plot one accuracy-runtime panel for each Report 2 method."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path


METHOD_ORDER = ("weno2", "weno3", "ader3")
BACKEND_ORDER = ("serial", "openmp", "cuda")
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
METHOD_LABELS = {
    "weno2": "WENO-2",
    "weno3": "WENO-3",
    "ader3": "ADER-3",
}
FIGURE_SIZE = (15.2, 5.6)
FIGURE_TITLE_SIZE = 21
PANEL_TITLE_SIZE = 18
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 15
LEGEND_SIZE = 15
ANNOTATION_SIZE = 11


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


def plain_log_axes(axis) -> None:
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    formatter = FuncFormatter(log_tick)
    for axis_obj in (axis.xaxis, axis.yaxis):
        axis_obj.set_major_locator(LogLocator(base=10))
        axis_obj.set_major_formatter(formatter)
        axis_obj.set_minor_formatter(NullFormatter())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument(
        "--time-metric",
        default="best_seconds",
        choices=("best_seconds", "mean_seconds", "median_seconds"),
    )
    parser.add_argument(
        "--error-metric",
        default="l1_rho",
        choices=("l1_rho", "l2_rho", "linf_rho"),
    )
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    time_metric_label = args.time_metric.replace("_seconds", "")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["text.parse_math"] = False
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    with args.csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit("accuracy-runtime CSV is empty")

    grouped: dict[tuple[str, str], list[tuple[int, float, float]]] = defaultdict(list)
    for row in rows:
        grouped[(row["method"], row["backend"])].append(
            (
                int(row["nx"]),
                float(row[args.time_metric]),
                float(row[args.error_metric]),
            )
        )
    for values in grouped.values():
        values.sort()

    methods = [method for method in METHOD_ORDER if any(row["method"] == method for row in rows)]
    fig, axes = plt.subplots(
        1,
        len(methods),
        figsize=FIGURE_SIZE,
        dpi=args.dpi,
        sharey=True,
    )
    if len(methods) == 1:
        axes = [axes]

    for axis, method in zip(axes, methods):
        for backend in BACKEND_ORDER:
            points = grouped.get((method, backend), [])
            if not points:
                continue
            axis.plot(
                [point[1] for point in points],
                [point[2] for point in points],
                color=COLORS[backend],
                marker=MARKERS[backend],
                linewidth=2,
                markersize=6,
            )
            if backend == "cuda":
                for resolution, runtime, error in points:
                    axis.annotate(
                        f"n={resolution}",
                        (runtime, error),
                        xytext=(5, 4),
                        textcoords="offset points",
                        fontsize=ANNOTATION_SIZE,
                        color="#333333",
                    )
        axis.set_title(
            METHOD_LABELS.get(method, method),
            fontsize=PANEL_TITLE_SIZE,
            pad=10,
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        plain_log_axes(axis)
        axis.set_xlabel(
            f"{time_metric_label} runtime (s)",
            fontsize=AXIS_LABEL_SIZE,
            labelpad=7,
        )
        axis.tick_params(axis="both", labelsize=TICK_LABEL_SIZE, pad=4)
        for spine in axis.spines.values():
            spine.set_linewidth(1.1)
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel(
        "L1 density error",
        fontsize=AXIS_LABEL_SIZE,
        labelpad=9,
    )

    handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[backend],
            marker=MARKERS[backend],
            linewidth=2,
            label=backend,
        )
        for backend in BACKEND_ORDER
    ]
    fig.suptitle(
        "Report 2: 2D Gaussian accuracy-runtime trade-off",
        y=0.985,
        fontsize=FIGURE_TITLE_SIZE,
    )
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.90),
        ncol=3,
        frameon=True,
        fontsize=LEGEND_SIZE,
    )
    fig.subplots_adjust(
        left=0.075, right=0.985, bottom=0.15, top=0.72, wspace=0.08
    )

    output_path = args.output or args.csv_path.with_suffix(
        f".{args.time_metric}.{args.error_metric}.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
