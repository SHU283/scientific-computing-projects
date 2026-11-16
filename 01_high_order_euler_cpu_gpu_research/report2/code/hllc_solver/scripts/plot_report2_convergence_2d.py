#!/usr/bin/env python3
"""Plot Report 2 smooth 2D convergence from one or more summary CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
LABELS = {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}
COLORS = {"weno2": "#d97706", "weno3": "#0f766e", "ader3": "#9333ea"}
MARKERS = {"weno2": "s", "weno3": "^", "ader3": "D"}
METRIC_LABELS = {
    "l1_rho": "L1 density error",
    "l2_rho": "L2 density error",
    "linf_rho": "Linf density error",
}


def log_tick(value: float, _position: int) -> str:
    if value <= 0:
        return ""
    exponent = round(math.log10(value))
    if not math.isclose(value, 10**exponent, rel_tol=1e-8, abs_tol=0.0):
        return ""
    if exponent == 0:
        return "1"
    return f"1e{exponent}"


def plain_log_y(axis) -> None:
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    axis.yaxis.set_major_locator(LogLocator(base=10))
    axis.yaxis.set_major_formatter(FuncFormatter(log_tick))
    axis.yaxis.set_minor_formatter(NullFormatter())


def load_rows(paths: list[Path]) -> dict[str, list[dict[str, float]]]:
    by_method: dict[str, dict[int, dict[str, float]]] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                method = row["method"].lower()
                if method not in METHODS:
                    continue
                n = int(row["n"])
                by_method.setdefault(method, {})[n] = {
                    "n": float(n),
                    "l1_rho": float(row["l1_rho"]),
                    "l2_rho": float(row["l2_rho"]),
                    "linf_rho": float(row["linf_rho"]),
                }
    return {
        method: [rows[n] for n in sorted(rows)]
        for method, rows in by_method.items()
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_paths",
        type=Path,
        nargs="*",
        help="Convergence summaries to merge by method and resolution.",
    )
    parser.add_argument(
        "--metric",
        choices=tuple(METRIC_LABELS),
        default="l1_rho",
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    paths = args.csv_paths or [
        root
        / "results"
        / "report2"
        / "validation_2d"
        / "report2_high_order_convergence_2d.csv",
        root
        / "results"
        / "report2"
        / "validation_2d_ader3"
        / "report2_high_order_convergence_2d.csv",
    ]
    paths = [path if path.is_absolute() else root / path for path in paths]
    for path in paths:
        if not path.exists():
            raise SystemExit(f"convergence summary not found: {path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["text.parse_math"] = False
    except ImportError:
        raise SystemExit("matplotlib is required")

    rows = load_rows(paths)
    if not rows:
        raise SystemExit("no Report 2 convergence rows found")

    fig, ax = plt.subplots(figsize=(8.2, 5.8), dpi=args.dpi)
    xticks: set[int] = set()
    for method in METHODS:
        method_rows = rows.get(method, [])
        if not method_rows:
            continue
        n_values = [int(row["n"]) for row in method_rows]
        errors = [row[args.metric] for row in method_rows]
        xticks.update(n_values)
        ax.plot(
            n_values,
            errors,
            color=COLORS[method],
            marker=MARKERS[method],
            linewidth=2.0,
            markersize=6.5,
            label=LABELS[method],
        )

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    plain_log_y(ax)
    ax.set_xticks(sorted(xticks))
    ax.set_xticklabels([str(value) for value in sorted(xticks)])
    ax.set_xlabel("Resolution n (n x n)")
    ax.set_ylabel(METRIC_LABELS[args.metric])
    ax.set_title("2D Gaussian-advection convergence")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output = args.output
    if output is None:
        output = (
            root
            / "results"
            / "report2"
            / "validation_2d"
            / f"report2_high_order_convergence_2d.{args.metric}.png"
        )
    elif not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
