#!/usr/bin/env python3
"""Plot Report 2 discontinuous runtime and speedup results."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import ticker

plt.rcParams["text.parse_math"] = False



METHODS = ("weno2", "weno3", "ader3")
METHOD_LABELS = {
    "weno2": "WENO-2",
    "weno3": "WENO-3",
    "ader3": "ADER-3",
}
BACKENDS = ("serial", "openmp", "cuda")
BACKEND_LABELS = {
    "serial": "serial",
    "openmp": "OpenMP",
    "cuda": "CUDA",
}
COLORS = {
    "serial": "#245783",
    "openmp": "#dd7500",
    "cuda": "#14827a",
}
MARKERS = {
    "serial": "o",
    "openmp": "s",
    "cuda": "^",
}
CASE_LABELS = {
    "quadrant_riemann_2d": "Quadrant Riemann",
    "shock_bubble_2d": "Shock-bubble",
}
FIGURE_TITLE_SIZE = 21
PANEL_TITLE_SIZE = 18
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 15
LEGEND_SIZE = 15


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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def grouped_rows(rows: list[dict[str, str]]):
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["problem"], row["method"], row["backend"])].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: int(row["nx"]))
    return grouped


def runtime_plot(
    rows: list[dict[str, str]],
    output: Path,
    metric: str,
    metric_label: str,
) -> None:
    cases = [case for case in CASE_LABELS if any(row["problem"] == case for row in rows)]
    grouped = grouped_rows(rows)
    figure, axes = plt.subplots(
        len(cases),
        len(METHODS),
        figsize=(15.2, 4.7 * len(cases)),
        squeeze=False,
        sharex=True,
    )

    for row_index, case in enumerate(cases):
        for column_index, method in enumerate(METHODS):
            axis = axes[row_index][column_index]
            for backend in BACKENDS:
                values = grouped.get((case, method, backend), [])
                if not values:
                    continue
                axis.plot(
                    [int(row["nx"]) for row in values],
                    [float(row[metric]) for row in values],
                    color=COLORS[backend],
                    marker=MARKERS[backend],
                    linewidth=2,
                    markersize=6,
                    label=BACKEND_LABELS[backend],
                )
            axis.set_yscale("log")
            plain_log_y(axis)
            axis.grid(True, alpha=0.25)
            axis.set_title(
                METHOD_LABELS[method],
                fontsize=PANEL_TITLE_SIZE,
                pad=9,
            )
            axis.set_xlabel(
                "resolution n",
                fontsize=AXIS_LABEL_SIZE,
                labelpad=6,
            )
            axis.tick_params(axis="both", labelsize=TICK_LABEL_SIZE, pad=3)
            for spine in axis.spines.values():
                spine.set_linewidth(1.1)
            if column_index == 0:
                axis.set_ylabel(
                    f"{CASE_LABELS[case]}\n{metric_label} runtime (s)",
                    fontsize=AXIS_LABEL_SIZE,
                    labelpad=8,
                )

    handles, labels = axes[0][0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=len(labels),
        frameon=True,
        fontsize=LEGEND_SIZE,
    )
    figure.suptitle(
        "Report 2 discontinuous 2D runtime scaling",
        y=0.995,
        fontsize=FIGURE_TITLE_SIZE,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.89), h_pad=2.1, w_pad=1.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def speedup_plot(
    rows: list[dict[str, str]],
    output: Path,
    metric: str,
) -> None:
    cases = [case for case in CASE_LABELS if any(row["problem"] == case for row in rows)]
    grouped = grouped_rows(rows)
    figure, axes = plt.subplots(
        len(cases),
        len(METHODS),
        figsize=(15.2, 4.7 * len(cases)),
        squeeze=False,
        sharex=True,
    )

    for row_index, case in enumerate(cases):
        for column_index, method in enumerate(METHODS):
            axis = axes[row_index][column_index]
            serial_values = {
                int(row["nx"]): float(row[metric])
                for row in grouped.get((case, method, "serial"), [])
            }
            resolutions = sorted(serial_values)
            axis.axhline(1.0, color="#666666", linestyle="--", linewidth=1.2)
            for backend in ("openmp", "cuda"):
                values = {
                    int(row["nx"]): float(row[metric])
                    for row in grouped.get((case, method, backend), [])
                }
                common = [resolution for resolution in resolutions if resolution in values]
                if not common:
                    continue
                axis.plot(
                    common,
                    [serial_values[resolution] / values[resolution] for resolution in common],
                    color=COLORS[backend],
                    marker=MARKERS[backend],
                    linewidth=2,
                    markersize=6,
                    label=BACKEND_LABELS[backend],
                )
            axis.grid(True, alpha=0.25)
            axis.set_ylim(bottom=0.0)
            axis.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
            axis.set_title(
                METHOD_LABELS[method],
                fontsize=PANEL_TITLE_SIZE,
                pad=9,
            )
            axis.set_xlabel(
                "resolution n",
                fontsize=AXIS_LABEL_SIZE,
                labelpad=6,
            )
            axis.tick_params(axis="both", labelsize=TICK_LABEL_SIZE, pad=3)
            for spine in axis.spines.values():
                spine.set_linewidth(1.1)
            if column_index == 0:
                axis.set_ylabel(
                    f"{CASE_LABELS[case]}\nspeedup vs serial",
                    fontsize=AXIS_LABEL_SIZE,
                    labelpad=8,
                )

    handles, labels = axes[0][0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=len(labels),
        frameon=True,
        fontsize=LEGEND_SIZE,
    )
    figure.suptitle(
        "Report 2 discontinuous 2D backend speedup",
        y=0.995,
        fontsize=FIGURE_TITLE_SIZE,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.89), h_pad=2.1, w_pad=1.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path")
    parser.add_argument(
        "--metric",
        choices=("best_seconds", "median_seconds"),
        default="best_seconds",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    rows = read_rows(csv_path)
    stem = csv_path.with_suffix("")
    metric_name = "best" if args.metric == "best_seconds" else "median"
    runtime_path = stem.with_name(f"{stem.name}.{metric_name}.runtime.png")
    speedup_path = stem.with_name(f"{stem.name}.{metric_name}.speedup.png")
    runtime_plot(rows, runtime_path, args.metric, metric_name)
    speedup_plot(rows, speedup_path, args.metric)
    print(f"saved plot: {runtime_path}")
    print(f"saved plot: {speedup_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
