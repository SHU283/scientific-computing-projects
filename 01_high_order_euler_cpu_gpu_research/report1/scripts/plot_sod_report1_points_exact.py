#!/usr/bin/env python3
"""Plot Report 1 Sod profiles as discrete numerical points with exact reference."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


GAMMA = 1.4


def load_series(path: Path, field: str) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row[field]))
    return xs, ys


def pressure_function(p: float, rho_i: float, p_i: float) -> tuple[float, float]:
    a_i = math.sqrt(GAMMA * p_i / rho_i)
    if p > p_i:
        a = 2.0 / ((GAMMA + 1.0) * rho_i)
        b = (GAMMA - 1.0) / (GAMMA + 1.0) * p_i
        root = math.sqrt(a / (p + b))
        f = (p - p_i) * root
        fd = root * (1.0 - 0.5 * (p - p_i) / (p + b))
    else:
        exponent = (GAMMA - 1.0) / (2.0 * GAMMA)
        ratio = p / p_i
        f = 2.0 * a_i / (GAMMA - 1.0) * (ratio**exponent - 1.0)
        fd = (1.0 / (rho_i * a_i)) * ratio ** (-(GAMMA + 1.0) / (2.0 * GAMMA))
    return f, fd


def solve_star_state() -> tuple[float, float]:
    rho_l, u_l, p_l = 1.0, 0.0, 1.0
    rho_r, u_r, p_r = 0.125, 0.0, 0.1
    p = 0.5 * (p_l + p_r)
    for _ in range(30):
        f_l, fd_l = pressure_function(p, rho_l, p_l)
        f_r, fd_r = pressure_function(p, rho_r, p_r)
        p_new = p - (f_l + f_r + u_r - u_l) / (fd_l + fd_r)
        if p_new <= 0.0:
            p_new = 0.5 * p
        if abs(p_new - p) / (0.5 * (p_new + p)) < 1e-12:
            p = p_new
            break
        p = p_new
    f_l, _ = pressure_function(p, rho_l, p_l)
    f_r, _ = pressure_function(p, rho_r, p_r)
    u = 0.5 * (u_l + u_r + f_r - f_l)
    return p, u


def sod_exact_rho(x: float, t: float, x0: float = 0.5) -> float:
    rho_l, u_l, p_l = 1.0, 0.0, 1.0
    rho_r, u_r, p_r = 0.125, 0.0, 0.1
    a_l = math.sqrt(GAMMA * p_l / rho_l)
    a_r = math.sqrt(GAMMA * p_r / rho_r)
    p_star, u_star = solve_star_state()
    s = (x - x0) / t

    if s <= u_star:
        if p_star > p_l:
            shock_speed = u_l - a_l * math.sqrt(
                (GAMMA + 1.0) / (2.0 * GAMMA) * p_star / p_l
                + (GAMMA - 1.0) / (2.0 * GAMMA)
            )
            if s <= shock_speed:
                return rho_l
            return rho_l * (
                (p_star / p_l + (GAMMA - 1.0) / (GAMMA + 1.0))
                / ((GAMMA - 1.0) / (GAMMA + 1.0) * p_star / p_l + 1.0)
            )

        head = u_l - a_l
        a_star = a_l * (p_star / p_l) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
        tail = u_star - a_star
        if s <= head:
            return rho_l
        if s >= tail:
            return rho_l * (p_star / p_l) ** (1.0 / GAMMA)
        a = 2.0 / (GAMMA + 1.0) * (a_l + 0.5 * (GAMMA - 1.0) * (u_l - s))
        return rho_l * (a / a_l) ** (2.0 / (GAMMA - 1.0))

    if p_star > p_r:
        shock_speed = u_r + a_r * math.sqrt(
            (GAMMA + 1.0) / (2.0 * GAMMA) * p_star / p_r
            + (GAMMA - 1.0) / (2.0 * GAMMA)
        )
        if s >= shock_speed:
            return rho_r
        return rho_r * (
            (p_star / p_r + (GAMMA - 1.0) / (GAMMA + 1.0))
            / ((GAMMA - 1.0) / (GAMMA + 1.0) * p_star / p_r + 1.0)
        )

    head = u_r + a_r
    a_star = a_r * (p_star / p_r) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
    tail = u_star + a_star
    if s >= head:
        return rho_r
    if s <= tail:
        return rho_r * (p_star / p_r) ** (1.0 / GAMMA)
    a = 2.0 / (GAMMA + 1.0) * (a_r - 0.5 * (GAMMA - 1.0) * (u_r - s))
    return rho_r * (a / a_r) ** (2.0 / (GAMMA - 1.0))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_paths", type=Path, nargs=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--dpi", type=int, default=260)
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=args.dpi)
    colors = ["#1f4e79", "#d97706", "#0f766e"]
    markers = ["o", "s", "^"]
    labels = ["Godunov", "MUSCL", "WENO-2"]

    for idx, path in enumerate(args.csv_paths):
        x, rho = load_series(path, "rho")
        ax.scatter(
            x[::4],
            rho[::4],
            s=18,
            color=colors[idx],
            marker=markers[idx],
            label=labels[idx],
            alpha=0.9,
        )

    xs = [i / 1200.0 for i in range(1201)]
    exact = [sod_exact_rho(x, args.final_time) for x in xs]
    ax.plot(xs, exact, "k--", linewidth=1.8, label="exact Sod solution")

    ax.set_xlabel("x", fontsize=14)
    ax.set_ylabel("rho", fontsize=14)
    ax.set_title("Sod shock tube density at t = 0.2", fontsize=16)
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=11)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    print(f"saved plot: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

