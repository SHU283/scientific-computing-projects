#!/usr/bin/env python3
"""Plot Report 2 1D convergence summaries."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
LABELS = {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}
COLORS = {"weno2": "#d97706", "weno3": "#0f766e", "ader3": "#9333ea"}
MARKERS = {"weno2": "s", "weno3": "^", "ader3": "D"}
PROBLEM_LABELS = {
    "smooth_advection": "1D sinusoidal advection",
    "gaussian_advection": "1D Gaussian advection",
}
METRIC_LABELS = {
    "l1_rho": "L1 density error",
    "l2_rho": "L2 density error",
    "linf_rho": "Linf density error",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_path",
        type=Path,
        nargs="?",
        default=Path("results/report2/validation_1d/report2_high_order_convergence_1d.csv"),
    )
    parser.add_argument("--metric", choices=tuple(METRIC_LABELS), default="l1_rho")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    csv_path = args.csv_path if args.csv_path.is_absolute() else root / args.csv_path
    if not csv_path.exists():
        raise SystemExit(f"convergence summary not found: {csv_path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["axes.labelsize"] = 16
        plt.rcParams["axes.titlesize"] = 17
        plt.rcParams["xtick.labelsize"] = 14
        plt.rcParams["ytick.labelsize"] = 14
        plt.rcParams["legend.fontsize"] = 14
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    def sci_tick(value: float, _position: int) -> str:
        if value <= 0.0:
            return ""
        exponent = int(round(__import__("math").log10(value)))
        if abs(value / (10.0**exponent) - 1.0) < 1.0e-8:
            return f"1e{exponent}"
        return ""

    grouped: dict[str, dict[str, list[tuple[int, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with csv_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            problem = row["problem"]
            method = row["method"]
            if method not in METHODS:
                continue
            grouped[problem][method].append((int(row["nx"]), float(row[args.metric])))

    problems = [problem for problem in PROBLEM_LABELS if problem in grouped]
    if not problems:
        raise SystemExit("no Report 2 1D convergence rows found")

    fig, axes = plt.subplots(
        1,
        len(problems),
        figsize=(7.6 * len(problems), 5.9),
        dpi=args.dpi,
        sharey=True,
    )
    if len(problems) == 1:
        axes = [axes]

    for axis, problem in zip(axes, problems):
        xticks: set[int] = set()
        for method in METHODS:
            points = sorted(grouped[problem].get(method, []))
            if not points:
                continue
            xs = [point[0] for point in points]
            ys = [point[1] for point in points]
            xticks.update(xs)
            axis.plot(
                xs,
                ys,
                color=COLORS[method],
                marker=MARKERS[method],
                linewidth=2.6,
                markersize=7.6,
                label=LABELS[method],
            )
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.yaxis.set_major_formatter(FuncFormatter(sci_tick))
        axis.set_xticks(sorted(xticks))
        axis.set_xticklabels([str(value) for value in sorted(xticks)])
        axis.set_xlabel("Resolution n", fontsize=16)
        axis.set_title(PROBLEM_LABELS[problem], fontsize=17, pad=12)
        axis.tick_params(axis="both", labelsize=14)
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel(METRIC_LABELS[args.metric], fontsize=16)
    axes[-1].legend(loc="best", fontsize=14)
    fig.suptitle("Report 2: 1D high-order convergence", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.91))

    output = args.output
    if output is None:
        output = csv_path.with_suffix(f".{args.metric}.png")
    elif not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
