#!/usr/bin/env python3
"""
Compare multiple 1D CSV solution curves on one figure.

For report figures, numerical finite-volume values can be shown with
markers while an analytic smooth-advection reference is overlaid as a line.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_series(path: Path, field_name: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row[field_name]))
    return xs, ys


def guess_problem_from_name(stem: str) -> str | None:
    if "smooth_advection" in stem:
        return "smooth_advection"
    if "gaussian_advection" in stem:
        return "gaussian_advection"
    return None


def exact_rho(problem: str, x: float, final_time: float, a: float, x_min: float, x_max: float) -> float:
    length = x_max - x_min
    shifted = x - a * final_time
    while shifted < x_min:
        shifted += length
    while shifted >= x_max:
        shifted -= length

    if problem == "smooth_advection":
        phase_shift = 0.173 * length
        return 1.0 + 0.2 * math.sin(2.0 * math.pi * (shifted - x_min - phase_shift) / length)
    if problem == "gaussian_advection":
        center = x_min + 0.317 * length
        sigma = 0.08 * length
        # Use the shortest periodic distance so the reference remains correct
        # after the translated Gaussian wraps around the domain.
        delta = shifted - center
        if delta > 0.5 * length:
            delta -= length
        elif delta < -0.5 * length:
            delta += length
        r = delta / sigma
        return 1.0 + 0.2 * math.exp(-0.5 * r * r)
    raise ValueError(f"unsupported exact problem: {problem}")


def default_stage_d_output(csv_paths: list[Path], field: str) -> Path:
    stem = "__vs__".join(path.stem for path in csv_paths)
    return csv_paths[0].resolve().parent.parent / "stageD" / f"{stem}.{field}.png"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", type=Path, nargs="+")
    parser.add_argument("--field", default="rho", choices=["rho", "u", "p", "energy"])
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(8.0, 5.0))
    parser.add_argument("--with-exact", action="store_true")
    parser.add_argument("--final-time", type=float, default=0.3)
    parser.add_argument("--adv-vel", type=float, default=0.73)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--markers", action="store_true", help="plot numerical values as discrete points")
    parser.add_argument("--marker-every", type=int, default=8)
    parser.add_argument("--line-width", type=float, default=1.6)
    parser.add_argument("--marker-size", type=float, default=4.5)
    parser.add_argument("--title", default=None)
    parser.add_argument("--title-fontsize", type=int, default=16)
    parser.add_argument("--label-fontsize", type=int, default=14)
    parser.add_argument("--tick-fontsize", type=int, default=12)
    parser.add_argument("--legend-fontsize", type=int, default=11)
    args = parser.parse_args()

    for path in args.csv_paths:
        if not path.exists():
            raise SystemExit(f"file not found: {path}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)
    colors = ["#1f4e79", "#d97706", "#0f766e", "#b91c1c", "#6d28d9"]
    marker_cycle = ["o", "s", "^", "D", "P"]

    problem = None
    for path in args.csv_paths:
        problem = problem or guess_problem_from_name(path.stem)

    for idx, path in enumerate(args.csv_paths):
        x, y = load_series(path, args.field)
        color = colors[idx % len(colors)]
        if args.markers:
            ax.plot(
                x,
                y,
                linestyle="-",
                linewidth=args.line_width,
                marker=marker_cycle[idx % len(marker_cycle)],
                markersize=args.marker_size,
                markevery=max(args.marker_every, 1),
                color=color,
                label=path.stem,
            )
        else:
            ax.plot(x, y, linewidth=args.line_width, color=color, label=path.stem)

    if args.with_exact and args.field == "rho" and problem is not None:
        # Draw the analytic reference on a fine plotting grid. The numerical
        # curves remain discrete, but the exact solution is a continuous target.
        n_exact = 1000
        xs = [
            args.x_min + (args.x_max - args.x_min) * k / (n_exact - 1)
            for k in range(n_exact)
        ]
        exact = [
            exact_rho(problem, xi, args.final_time, args.adv_vel, args.x_min, args.x_max)
            for xi in xs
        ]
        ax.plot(xs, exact, "k--", linewidth=1.8, label=f"{problem} exact")

    ax.set_xlabel("x", fontsize=args.label_fontsize)
    ax.set_ylabel(args.field, fontsize=args.label_fontsize)
    ax.set_title(args.title or f"1D comparison [{args.field}]", fontsize=args.title_fontsize)
    ax.tick_params(axis="both", labelsize=args.tick_fontsize)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=args.legend_fontsize)
    fig.tight_layout()

    output_path = args.output or default_stage_d_output(args.csv_paths, args.field)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved comparison plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
