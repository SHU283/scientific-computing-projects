#!/usr/bin/env python3
"""
Plot 1D field from solver CSV, optionally comparing with exact solution.
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
        r = (shifted - center) / sigma
        return 1.0 + 0.2 * math.exp(-0.5 * r * r)
    raise ValueError(f"unsupported exact problem: {problem}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--field", default="rho", choices=["rho", "u", "p", "energy"])
    parser.add_argument("--final-time", type=float, default=0.3)
    parser.add_argument("--adv-vel", type=float, default=0.73)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--with-exact", action="store_true")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=240)
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise SystemExit(f"file not found: {args.csv_path}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required. Try: pip install matplotlib")

    x, y = load_series(args.csv_path, args.field)
    problem = guess_problem_from_name(args.csv_path.stem)

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=args.dpi)
    ax.plot(x, y, linewidth=2.0, label=args.csv_path.stem)

    if args.with_exact and args.field == "rho" and problem is not None:
        exact = [
            exact_rho(problem, xi, args.final_time, args.adv_vel, args.x_min, args.x_max)
            for xi in x
        ]
        ax.plot(x, exact, "--", linewidth=1.6, label=f"{problem} exact")

    ax.set_xlabel("x")
    ax.set_ylabel(args.field)
    ax.set_title(f"{args.csv_path.name} [{args.field}]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    if args.output is None:
        output_path = args.csv_path.with_suffix(f".{args.field}.png")
    else:
        output_path = args.output

    fig.savefig(output_path, bbox_inches="tight")
    print(f"saved plot: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
