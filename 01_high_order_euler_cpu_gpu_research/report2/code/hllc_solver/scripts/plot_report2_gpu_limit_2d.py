#!/usr/bin/env python3
"""Plot the Report 2 CUDA-only high-resolution smooth sweep."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path


METHOD_ORDER = ("weno2", "weno3", "ader3")
# Draw ADER-3 before WENO-3 because their density errors are very close on the
# high-resolution Gaussian sweep; otherwise WENO-3 can disappear underneath it.
PLOT_ORDER = ("weno2", "ader3", "weno3")
METHOD_LABELS = {
    "weno2": "WENO-2",
    "weno3": "WENO-3",
    "ader3": "ADER-3",
}
COLORS = {
    "weno2": "#d97706",
    "weno3": "#0f766e",
    "ader3": "#7c3aed",
}
MARKERS = {
    "weno2": "s",
    "weno3": "^",
    "ader3": "D",
}
LINESTYLES = {
    "weno2": "-",
    "weno3": "-",
    "ader3": "-",
}


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
    return f"1e{exponent}"


def apply_log_format(axis, *, x: bool = False, y: bool = False) -> None:
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    if x:
        axis.xaxis.set_major_locator(LogLocator(base=10))
        axis.xaxis.set_major_formatter(FuncFormatter(log_tick))
        axis.xaxis.set_minor_formatter(NullFormatter())
    if y:
        axis.yaxis.set_major_locator(LogLocator(base=10))
        axis.yaxis.set_major_formatter(FuncFormatter(log_tick))
        axis.yaxis.set_minor_formatter(NullFormatter())


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"empty CSV: {path}")
    return rows


def grouped_points(
    rows: list[dict[str, str]],
    x_field: str,
    y_field: str,
) -> dict[str, list[tuple[float, float]]]:
    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in rows:
        if row["backend"] != "cuda":
            continue
        grouped[row["method"]].append((float(row[x_field]), float(row[y_field])))
    for values in grouped.values():
        values.sort()
    return grouped


def resolution_ticks(resolutions: list[int]) -> list[int]:
    """Keep the physical resolution axis readable on long GPU-limit sweeps."""
    if resolutions[-1] >= 7000 and resolutions[0] < 1000:
        return list(range(1000, 7001, 1000))
    if len(resolutions) <= 8:
        return resolutions
    indices = (0, 3, 6, 9, len(resolutions) - 1)
    return [resolutions[index] for index in indices]


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
    args = parser.parse_args()

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams["font.family"] = "Arial"
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    rows = load_rows(args.csv_path)
    resolutions = sorted({int(row["nx"]) for row in rows if row["backend"] == "cuda"})
    display_resolutions = resolution_ticks(resolutions)
    runtime = grouped_points(rows, "nx", args.time_metric)
    error = grouped_points(rows, "nx", args.error_metric)
    tradeoff = grouped_points(rows, args.time_metric, args.error_metric)

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.2), dpi=args.dpi)
    panels = [
        (axes[0], runtime, "resolution n", "CUDA runtime (s)", "Runtime scaling", "linear", "log"),
        (axes[1], error, "resolution n", "L1 density error", "Error refinement", "linear", "log"),
        (
            axes[2],
            tradeoff,
            "CUDA runtime (s)",
            "L1 density error",
            "Accuracy-runtime path",
            "log",
            "log",
        ),
    ]

    for axis, grouped, xlabel, ylabel, title, xscale, yscale in panels:
        for method in PLOT_ORDER:
            points = grouped.get(method, [])
            if not points:
                continue
            xs = [point[0] for point in points]
            ys = [point[1] for point in points]
            axis.plot(
                xs,
                ys,
                color=COLORS[method],
                marker=MARKERS[method],
                linestyle=LINESTYLES[method],
                linewidth=1.8,
                markersize=5.5,
                markeredgecolor="white",
                markeredgewidth=0.45,
                zorder=3 if method == "weno3" else 2,
                label=METHOD_LABELS[method],
            )
        axis.set_title(title, fontsize=18, pad=11)
        axis.set_xlabel(xlabel, fontsize=16, labelpad=8)
        axis.set_ylabel(ylabel, fontsize=16, labelpad=9)
        axis.tick_params(axis="both", labelsize=15)
        axis.set_xscale(xscale)
        axis.set_yscale(yscale)
        apply_log_format(axis, x=(xscale == "log"), y=(yscale == "log"))
        if xscale == "linear":
            axis.set_xticks(display_resolutions)
        axis.grid(True, which="both", alpha=0.25)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.93), ncol=3, frameon=True)
    fig.suptitle(
        "Report 2: CUDA-only high-resolution 2D Gaussian sweep",
        y=0.995,
        fontsize=21,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.84))

    output = args.csv_path.with_suffix(
        f".{args.time_metric}.{args.error_metric}.gpu_limit.png"
    )
    fig.savefig(output, bbox_inches="tight")
    print(f"saved plot: {output}")

    for grouped, suffix, title in [
        (runtime, f".{args.time_metric}.gpu_runtime.png", "Runtime scaling"),
        (error, f".{args.error_metric}.gpu_error.png", "Error refinement"),
        (
            tradeoff,
            f".{args.time_metric}.{args.error_metric}.gpu_tradeoff.png",
            "Accuracy-runtime path",
        ),
    ]:
        single, single_axis = plt.subplots(figsize=(8.8, 5.2), dpi=args.dpi)
        for method in PLOT_ORDER:
            points = grouped.get(method, [])
            if not points:
                continue
            single_axis.plot(
                [point[0] for point in points],
                [point[1] for point in points],
                color=COLORS[method],
                marker=MARKERS[method],
                linestyle=LINESTYLES[method],
                linewidth=1.8,
                markersize=5.5,
                markeredgecolor="white",
                markeredgewidth=0.45,
                zorder=3 if method == "weno3" else 2,
                label=METHOD_LABELS[method],
            )
        if "runtime" in suffix or "error" in suffix:
            single_axis.set_xlabel("resolution n")
            single_axis.set_xticks(display_resolutions)
        else:
            single_axis.set_xlabel("CUDA runtime (s)")
            single_axis.set_xscale("log")
        single_axis.set_ylabel(
            "CUDA runtime (s)" if "runtime" in suffix else "L1 density error",
            fontsize=16,
        )
        single_axis.tick_params(axis="both", labelsize=15)
        single_axis.set_title(title, fontsize=18, pad=11)
        single_axis.set_yscale("log")
        apply_log_format(single_axis, x=("tradeoff" in suffix), y=True)
        single_axis.grid(True, which="both", alpha=0.25)
        single_axis.legend(frameon=True)
        single.tight_layout()
        single_path = args.csv_path.with_suffix(suffix)
        single.savefig(single_path, bbox_inches="tight")
        print(f"saved plot: {single_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
