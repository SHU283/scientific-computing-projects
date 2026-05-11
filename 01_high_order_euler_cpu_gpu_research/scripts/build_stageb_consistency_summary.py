#!/usr/bin/env python3
"""
Summarize Stage B backend consistency checks for the Report 1 method set.

Outputs:
1. CSV with per-field and overall max absolute differences vs serial
2. Clean summary figure grouped by comparison backend (openmp / cuda)
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


METHODS = ["godunov", "muscl", "weno2"]
PROBLEMS = [
    "smooth_advection",
    "gaussian_advection",
    "sod",
    "gaussian_advection_2d",
    "quadrant_riemann_2d",
    "shock_bubble_2d",
]
BACKENDS = ["openmp", "cuda"]


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def row_key(row: dict[str, str]) -> tuple[str, ...]:
    if "x" in row:
        return (row["x"],)
    if "i" in row and "j" in row:
        return (row["i"], row["j"])
    raise ValueError("unsupported csv coordinate format")


def compare_csvs(serial_path: Path, other_path: Path) -> dict[str, float]:
    rows_a = load_rows(serial_path)
    rows_b = load_rows(other_path)
    if not rows_a:
        raise ValueError(f"empty csv: {serial_path.name}")

    map_a = {row_key(row): row for row in rows_a}
    map_b = {row_key(row): row for row in rows_b}
    shared_keys = [key for key in map_a.keys() if key in map_b]
    if not shared_keys:
        raise ValueError(f"no shared coordinates: {serial_path.name} vs {other_path.name}")

    fields = [key for key in rows_a[0].keys() if key not in {"x", "i", "j"}]
    max_diff = {field: 0.0 for field in fields}

    for key in shared_keys:
        row_a = map_a[key]
        row_b = map_b[key]
        for field in fields:
            diff = abs(float(row_a[field]) - float(row_b[field]))
            if diff > max_diff[field]:
                max_diff[field] = diff

    max_diff["max_overall"] = max(max_diff.values()) if max_diff else 0.0
    return max_diff


def problem_label(problem: str) -> str:
    mapping = {
        "smooth_advection": "1D smooth",
        "gaussian_advection": "1D gaussian",
        "sod": "1D sod",
        "gaussian_advection_2d": "2D smooth",
        "quadrant_riemann_2d": "2D quadrant",
        "shock_bubble_2d": "2D shock-bubble",
    }
    return mapping.get(problem, problem.replace("_", " "))


def display_label(problem: str, method: str) -> str:
    return f"{problem_label(problem)} / {method}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stageb", default="results/stageB")
    parser.add_argument("--output-csv", default="results/stageB/backend_consistency_report1_summary.csv")
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(16.0, 9.0))
    parser.add_argument("--output-plot", default="results/stageB/backend_consistency_report1_summary.max_delta.png")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    stageb_dir = root / args.stageb
    csv_out = root / args.output_csv
    plot_out = root / args.output_plot

    rows: list[dict[str, str | float]] = []
    for problem in PROBLEMS:
        for method in METHODS:
            serial_path = stageb_dir / f"{problem}_serial_{method}.csv"
            if not serial_path.exists():
                continue
            for backend in BACKENDS:
                other_path = stageb_dir / f"{problem}_{backend}_{method}.csv"
                if not other_path.exists():
                    continue
                diffs = compare_csvs(serial_path, other_path)
                rows.append(
                    {
                        "problem": problem,
                        "method": method,
                        "backend": backend,
                        "rho": diffs.get("rho", 0.0),
                        "u": diffs.get("u", 0.0),
                        "v": diffs.get("v", 0.0),
                        "p": diffs.get("p", 0.0),
                        "energy": diffs.get("energy", 0.0),
                        "max_overall": diffs["max_overall"],
                    }
                )

    if not rows:
        raise SystemExit("no Stage B comparison pairs were found")

    csv_out.parent.mkdir(parents=True, exist_ok=True)
    with csv_out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["problem", "method", "backend", "rho", "u", "v", "p", "energy", "max_overall"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"saved summary csv: {csv_out}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.lines import Line2D
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        raise SystemExit("matplotlib and numpy are required. Try: pip install matplotlib numpy")

    plt.rcParams["text.parse_math"] = False

    backend_rows = {backend: [row for row in rows if row["backend"] == backend] for backend in BACKENDS}
    field_colors = {"u": "#2563eb", "v": "#7c3aed", "energy": "#d97706", "rho": "#059669", "p": "#dc2626"}

    def dominant_field(row: dict[str, str | float]) -> tuple[str, float]:
        candidates = ["rho", "u", "v", "p", "energy"]
        best_field = "rho"
        best_value = -1.0
        for field in candidates:
            value = float(row[field])
            if value > best_value:
                best_value = value
                best_field = field
        return best_field, best_value

    def format_sci(value: float, _pos: int) -> str:
        if value <= 0.0:
            return "0"
        return f"{value:.0e}".replace("e-0", "e-").replace("e+0", "e+")

    fig_width, fig_height = tuple(args.figsize)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(max(fig_width, 14.0), max(fig_height, 9.8)),
        dpi=args.dpi,
        constrained_layout=False,
        gridspec_kw={"height_ratios": [1.0, 2.2]},
    )

    openmp_ax, cuda_ax = axes
    openmp_rows = backend_rows["openmp"]
    openmp_nonzero = [row for row in openmp_rows if float(row["max_overall"]) > 0.0]

    openmp_ax.set_title("vs serial: openmp", fontsize=14)
    openmp_ax.axis("off")
    if not openmp_nonzero:
        openmp_ax.text(
            0.5,
            0.60,
            "All OpenMP comparisons exactly match serial",
            ha="center",
            va="center",
            fontsize=15,
            color="#b45309",
            weight="bold",
            transform=openmp_ax.transAxes,
        )
        openmp_ax.text(
            0.5,
            0.38,
            f"Checked {len(openmp_rows)} report-1 pairs; max abs difference = 0.0",
            ha="center",
            va="center",
            fontsize=12,
            color="#666666",
            transform=openmp_ax.transAxes,
        )
    else:
        openmp_ax.text(
            0.5,
            0.5,
            "Unexpected nonzero OpenMP differences detected",
            ha="center",
            va="center",
            fontsize=14,
            color="#b91c1c",
            transform=openmp_ax.transAxes,
        )

    cuda_nonzero = [row for row in backend_rows["cuda"] if float(row["max_overall"]) > 0.0]
    cuda_nonzero.sort(key=lambda row: float(row["max_overall"]), reverse=True)

    labels = [display_label(str(row["problem"]), str(row["method"])) for row in cuda_nonzero]
    values = [float(row["max_overall"]) for row in cuda_nonzero]
    y = np.arange(len(cuda_nonzero))

    xmin = 1e-16
    xmax = max(values) * 2.5 if values else 1e-5
    cuda_ax.set_title("vs serial: cuda", fontsize=14)
    cuda_ax.set_xscale("log", base=10)
    cuda_ax.set_xlim(xmin, xmax)
    cuda_ax.set_xlabel("max abs difference", fontsize=11)
    cuda_ax.set_ylabel("problem / method", fontsize=11)
    cuda_ax.set_yticks(y)
    cuda_ax.set_yticklabels(labels, fontsize=10)
    cuda_ax.tick_params(axis="x", labelsize=10)
    cuda_ax.tick_params(axis="y", labelsize=10)
    cuda_ax.grid(True, axis="x", which="major", alpha=0.25)
    cuda_ax.grid(True, axis="x", which="minor", alpha=0.08)
    cuda_ax.invert_yaxis()
    cuda_ax.xaxis.set_major_formatter(FuncFormatter(format_sci))

    for yi, row, value in zip(y, cuda_nonzero, values):
        field, _ = dominant_field(row)
        color = field_colors.get(field, "#333333")
        cuda_ax.hlines(yi, xmin, value, color=color, alpha=0.55, linewidth=2.2)
        cuda_ax.plot(value, yi, marker="o", color=color, markersize=8)
        cuda_ax.text(
            min(value * 1.14, xmax * 0.96),
            yi,
            f"{field}: {value:.1e}",
            va="center",
            ha="left",
            fontsize=9,
            color=color,
        )

    legend_handles = [
        Line2D([0], [0], marker="o", color=field_colors[field], linewidth=2.2, markersize=7, label=field)
        for field in ["u", "v", "energy"]
        if any(dominant_field(row)[0] == field for row in cuda_nonzero)
    ]
    if legend_handles:
        cuda_ax.legend(handles=legend_handles, title="dominant field", loc="lower right", fontsize=9, title_fontsize=10)

    fig.suptitle("Stage B Backend Consistency Summary (Report 1 methods)", fontsize=16)
    fig.subplots_adjust(left=0.26, right=0.97, bottom=0.08, top=0.92, hspace=0.22)
    fig.savefig(plot_out)
    print(f"saved summary plot: {plot_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
