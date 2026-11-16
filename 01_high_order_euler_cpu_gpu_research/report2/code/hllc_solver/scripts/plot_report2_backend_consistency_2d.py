#!/usr/bin/env python3
"""Plot Report 2 backend consistency differences against serial results."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FIELDS = ("rho", "u", "v", "p", "energy")
FIELD_COLORS = {
    "rho": "#1f4e79",
    "u": "#2563eb",
    "v": "#7c3aed",
    "p": "#d97706",
    "energy": "#0f766e",
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


def plain_log_x(axis) -> None:
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    axis.xaxis.set_major_locator(LogLocator(base=10))
    axis.xaxis.set_major_formatter(FuncFormatter(log_tick))
    axis.xaxis.set_minor_formatter(NullFormatter())


def label(problem: str, method: str) -> str:
    names = {
        "gaussian_advection_2d": "2D Gaussian",
        "quadrant_riemann_2d": "2D quadrant",
        "shock_bubble_2d": "2D shock-bubble",
    }
    return f"{names.get(problem, problem)} / {method.upper()}"


def dominant(row: dict[str, str]) -> tuple[str, float]:
    values = {
        field: float(row[f"max_abs_{field}_vs_serial"])
        for field in FIELDS
    }
    return max(values.items(), key=lambda item: item[1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv_path",
        type=Path,
        nargs="?",
        default=Path(
            "results/report2/backend_consistency_2d_lovelace/"
            "backend_consistency_2d_report2_lovelace_a30.csv"
        ),
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=260)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    csv_path = args.csv_path if args.csv_path.is_absolute() else root / args.csv_path
    if not csv_path.exists():
        raise SystemExit(f"backend consistency CSV not found: {csv_path}")

    rows: list[dict[str, str]] = []
    with csv_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["backend"] == "serial":
                continue
            field, value = dominant(row)
            row["dominant_field"] = field
            row["dominant_value"] = value
            rows.append(row)

    cuda_rows = [row for row in rows if row["backend"] == "cuda"]
    cuda_rows.sort(key=lambda row: float(row["dominant_value"]), reverse=True)
    openmp_max = max(
        (float(row["dominant_value"]) for row in rows if row["backend"] == "openmp"),
        default=0.0,
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        plt.rcParams["text.parse_math"] = False
    except ImportError as exc:
        raise SystemExit("matplotlib is required") from exc

    fig, ax = plt.subplots(figsize=(9.4, 5.8), dpi=args.dpi)
    xmin = 1e-15
    xmax = max(float(row["dominant_value"]) for row in cuda_rows) * 3.0
    y_positions = list(range(len(cuda_rows)))

    for y, row in zip(y_positions, cuda_rows):
        field = str(row["dominant_field"])
        value = float(row["dominant_value"])
        color = FIELD_COLORS[field]
        ax.hlines(y, xmin, value, color=color, linewidth=2.2, alpha=0.6)
        ax.plot(value, y, marker="o", markersize=7.5, color=color)
        text_x = min(value * 1.16, xmax * 0.93)
        ha = "left" if text_x < xmax * 0.90 else "right"
        ax.text(text_x, y, f"{field}: {value:.1e}", va="center", ha=ha, fontsize=9, color=color)

    ax.set_xscale("log")
    plain_log_x(ax)
    ax.set_xlim(xmin, xmax)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([label(row["problem"], row["method"]) for row in cuda_rows], fontsize=9)
    ax.invert_yaxis()
    ax.grid(True, axis="x", which="major", alpha=0.25)
    ax.grid(True, axis="x", which="minor", alpha=0.08)
    ax.set_xlabel("maximum absolute difference vs serial")
    ax.set_ylabel("problem / method")
    ax.set_title("Report 2 backend consistency on Lovelace A30")
    legend_fields = sorted({str(row["dominant_field"]) for row in cuda_rows})
    handles = [
        Line2D([0], [0], marker="o", color=FIELD_COLORS[field], linewidth=2.2, markersize=7, label=field)
        for field in legend_fields
    ]
    ax.legend(handles=handles, title="dominant field", loc="lower right")
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 1.0))
    fig.text(
        0.12,
        0.025,
        f"OpenMP maximum difference vs serial: {openmp_max:.1e} in exported fields.",
        ha="left",
        va="bottom",
        fontsize=10,
    )

    output = args.output
    if output is None:
        output = csv_path.with_suffix(".max_delta.png")
    output = output if output.is_absolute() else root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
