#!/usr/bin/env python3
"""
Plot accuracy vs runtime for the selected problem.
Originally added for 1D high-order tests, but it is also useful for
Report 1's 2D smooth-method comparison table.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def format_log_tick(value: float) -> str:
    if value <= 0.0:
        return ""
    if value >= 1.0:
        return f"{value:g}"
    if value >= 1e-2:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    if value >= 1e-3:
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return f"{value:.1e}"


def merge_label_styles(
    base: dict[tuple[str, str, str], dict[str, object]],
    override: dict[tuple[str, str, str], dict[str, object]],
) -> dict[tuple[str, str, str], dict[str, object]]:
    merged = dict(base)
    merged.update(override)
    return merged


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--problem", default=None)
    parser.add_argument("--error-metric", default="l1_rho", choices=["l1_rho", "l2_rho", "linf_rho"])
    parser.add_argument("--time-metric", default="best_seconds", choices=["best_seconds", "mean_seconds", "median_seconds"])
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 6.0))
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--no-point-labels", action="store_true")
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, NullFormatter
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    plt.rcParams["text.parse_math"] = False

    def pretty_label(text: str) -> str:
        return text.replace("_", " ")

    rows = load_rows(args.csv_path)
    if args.problem is not None:
        rows = [row for row in rows if row["problem"] == args.problem]
    if not rows:
        raise SystemExit("no rows matched the requested plot selection")

    colors = {"serial": "#1f4e79", "openmp": "#d97706", "cuda": "#0f766e"}
    markers = {
        "godunov": "o",
        "muscl": "s",
        "weno2": "^",
        "weno3": "D",
        "ader3": "P",
    }

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)

    label_styles = {
        ("godunov", "serial", "50"): {"xytext": (0, 6), "ha": "left", "va": "bottom"},
        ("godunov", "openmp", "50"): {"xytext": (4, -8), "ha": "left", "va": "top"},
        ("godunov", "cuda", "50"): {"xytext": (4, 4), "ha": "left", "va": "bottom"},
        ("godunov", "serial", "100"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
        ("godunov", "openmp", "100"): {"xytext": (6, 14), "ha": "left", "va": "bottom"},
        ("godunov", "cuda", "100"): {"xytext": (8, 4), "ha": "left", "va": "bottom"},
        ("muscl", "serial", "50"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
        ("muscl", "openmp", "50"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
        ("muscl", "cuda", "50"): {"xytext": (6, 6), "ha": "left", "va": "bottom"},
        ("muscl", "serial", "100"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
        ("muscl", "openmp", "100"): {"xytext": (6, 8), "ha": "left", "va": "bottom"},
        ("muscl", "cuda", "100"): {"xytext": (6, -8), "ha": "left", "va": "top"},
        ("weno2", "serial", "50"): {"xytext": (6, -8), "ha": "left", "va": "top"},
        ("weno2", "openmp", "50"): {"xytext": (6, -8), "ha": "left", "va": "top"},
        ("weno2", "cuda", "50"): {"xytext": (6, 6), "ha": "left", "va": "bottom"},
        ("weno2", "serial", "100"): {"xytext": (-6, 4), "ha": "right", "va": "bottom"},
        ("weno2", "openmp", "100"): {"xytext": (6, 8), "ha": "left", "va": "bottom"},
        ("weno2", "cuda", "100"): {"xytext": (6, -8), "ha": "left", "va": "top"},
    }
    if args.error_metric in {"l1_rho", "l2_rho"}:
        label_styles[("godunov", "cuda", "100")] = {"xytext": (8, -8), "ha": "left", "va": "top"}

    if args.problem == "gaussian_advection":
        label_styles = merge_label_styles(
            label_styles,
            {
                ("muscl", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("muscl", "openmp", "400"): {"xytext": (-6, 4), "ha": "right", "va": "bottom"},
                ("weno2", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("weno2", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
                ("weno3", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("weno3", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
                ("ader3", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("ader3", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
            },
        )
    elif args.problem == "smooth_advection":
        label_styles = merge_label_styles(
            label_styles,
            {
                ("muscl", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("muscl", "openmp", "400"): {"xytext": (-6, 4), "ha": "right", "va": "bottom"},
                ("weno2", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("weno2", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
                ("weno3", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("weno3", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
                ("ader3", "serial", "400"): {"xytext": (6, 4), "ha": "left", "va": "bottom"},
                ("ader3", "openmp", "400"): {"xytext": (6, -8), "ha": "left", "va": "top"},
            },
        )

    for row in rows:
        x = float(row[args.time_metric])
        y = float(row[args.error_metric])
        backend = row["backend"]
        method = row["method"]
        nx = row["nx"]

        ax.scatter(
            x,
            y,
            s=80,
            color=colors.get(backend, "#333333"),
            marker=markers.get(method, "o"),
            alpha=0.9,
        )
        if not args.no_point_labels:
            style = label_styles.get((method, backend, nx), {"xytext": (6, 4), "ha": "left", "va": "bottom"})
            ax.annotate(
                f"{backend}-{method}\nn={nx}",
                (x, y),
                textcoords="offset points",
                xytext=style["xytext"],
                fontsize=6.5,
                ha=style["ha"],
                va=style["va"],
            )

    title_problem = args.problem if args.problem is not None else "all"
    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)
    ax.set_xlabel(f"{pretty_label(args.time_metric)} (s)")
    ax.set_ylabel(pretty_label(args.error_metric))
    ax.set_title(f"Accuracy vs Runtime ({pretty_label(title_problem)})")
    ax.grid(True, alpha=0.25)
    ax.margins(x=0.06, y=0.14)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: format_log_tick(value)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: format_log_tick(value)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    from matplotlib.lines import Line2D

    present_backends = []
    for row in rows:
        backend = row["backend"]
        if backend not in present_backends:
            present_backends.append(backend)

    present_methods = []
    for row in rows:
        method = row["method"]
        if method not in present_methods:
            present_methods.append(method)

    backend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors.get(name, "#333333"), markersize=8, label=name)
        for name in present_backends
    ]
    method_handles = [
        Line2D([0], [0], marker=markers.get(name, "o"), color="#444444", linestyle="None", markersize=8, label=name)
        for name in present_methods
    ]
    legend1 = ax.legend(handles=backend_handles, title="backend", loc="upper left")
    ax.add_artist(legend1)
    ax.legend(handles=method_handles, title="method", loc="lower right")

    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.12, top=0.86)

    if args.output is None:
        suffix_problem = "all" if args.problem is None else args.problem
        output_path = args.csv_path.with_suffix(f".{suffix_problem}.{args.time_metric}.{args.error_metric}.png")
    else:
        output_path = args.output

    fig.savefig(output_path)
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
