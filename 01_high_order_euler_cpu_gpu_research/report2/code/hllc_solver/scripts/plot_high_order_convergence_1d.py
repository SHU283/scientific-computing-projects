#!/usr/bin/env python3
"""
Plot 1D high-order convergence curves for a chosen smooth problem.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


METHODS = ["godunov", "muscl", "weno2", "weno3", "ader3"]


def infer_title(problem: str, methods: list[str]) -> str:
    pretty_problem = problem.replace("_", " ")
    if methods == ["muscl", "weno2"]:
        return f"1D Convergence ({pretty_problem}, second-order methods)"
    if methods == ["godunov", "muscl", "weno2"]:
        return f"1D Convergence ({pretty_problem}, Report 1 methods)"
    if methods == METHODS:
        return f"1D High-Order Convergence ({pretty_problem})"
    method_text = ", ".join(methods)
    return f"1D Convergence ({pretty_problem}; {method_text})"


def load_rows(path: Path) -> list[dict[str, float]]:
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows: list[dict[str, float]] = []
        for row in reader:
            rows.append(
                {
                    "nx": float(row["nx"]),
                    "l1_rho": float(row["l1_rho"]),
                    "l2_rho": float(row["l2_rho"]),
                    "linf_rho": float(row["linf_rho"]),
                }
            )
    return rows


def metric_label(metric: str) -> str:
    mapping = {
        "l1_rho": "L1(rho) error",
        "l2_rho": "L2(rho) error",
        "linf_rho": "Linf(rho) error",
    }
    return mapping.get(metric, metric.replace("_", " "))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--problem", default="gaussian_advection", choices=["gaussian_advection", "smooth_advection"])
    parser.add_argument("--metric", default="l1_rho", choices=["l1_rho", "l2_rho", "linf_rho"])
    parser.add_argument("--methods", nargs="+", default=METHODS, choices=METHODS)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--title", default=None)
    parser.add_argument("--dpi", type=int, default=240)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent

    try:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    plt.rcParams["text.parse_math"] = False

    colors = {
        "godunov": "#6b7280",
        "muscl": "#1f4e79",
        "weno2": "#d97706",
        "weno3": "#0f766e",
        "ader3": "#9333ea",
    }
    markers = {
        "godunov": "o",
        "muscl": "s",
        "weno2": "^",
        "weno3": "o",
        "ader3": "D",
    }

    fig, ax = plt.subplots(figsize=(8.0, 6.0), dpi=args.dpi)

    plotted_any = False
    xtick_values: set[float] = set()
    for method in args.methods:
        path = root / "results" / "stageA" / f"{args.problem}_1d_{method}_convergence.csv"
        if not path.exists():
            continue
        rows = load_rows(path)
        if not rows:
            continue
        plotted_any = True
        x = [row["nx"] for row in rows]
        y = [row[args.metric] for row in rows]
        xtick_values.update(x)
        ax.plot(
            x,
            y,
            marker=markers[method],
            linewidth=2.0,
            markersize=6,
            color=colors[method],
            label=method,
        )

    if not plotted_any:
        raise SystemExit("no convergence csv files found for the requested problem")

    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=10)
    ax.set_xlabel("nx")
    ax.set_ylabel(metric_label(args.metric))
    ax.set_title(args.title or infer_title(args.problem, args.methods))
    ax.grid(True, alpha=0.25)
    ax.legend()
    xticks = sorted(xtick_values)
    if xticks:
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(int(x)) for x in xticks])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.90)

    if args.output is None:
        suffix = "high_order" if args.methods == METHODS else "__".join(args.methods)
        output_path = root / "results" / "stageA" / f"{args.problem}_1d_{suffix}.{args.metric}.png"
    else:
        output_path = args.output
    fig.savefig(output_path)
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
