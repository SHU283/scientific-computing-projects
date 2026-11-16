from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import ticker


ROOT = Path(__file__).resolve().parents[3]
CSV = (
    ROOT
    / "final_project"
    / "results"
    / "report2"
    / "performance"
    / "discontinuous_refresh_hll_20260705"
    / "runtime_discontinuous_2d_report2_hll_lovelace_a30_25_50_100_200_600_refresh.csv"
)
OUT = ROOT / "analysis_outputs" / "hll_vs_hllc" / "figures" / "hll_vs_hllc"

PROBLEM_LABELS = {
    "quadrant_riemann_2d": "Quadrant Riemann",
    "shock_bubble_2d": "Shock-bubble",
}
METHOD_LABELS = {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}
BACKEND_LABELS = {"serial": "serial", "openmp": "OpenMP", "cuda": "CUDA"}
COLORS = {"serial": "#1f4e79", "openmp": "#e07a00", "cuda": "#0f7f76"}
MARKERS = {"serial": "o", "openmp": "s", "cuda": "^"}
RESOLUTIONS = [25, 50, 100, 200, 600]


def plain_log_tick(value: float, _pos: int) -> str:
    if value <= 0:
        return ""
    if value >= 1000:
        return f"{value / 1000:g}k"
    if value >= 1:
        return f"{value:g}"
    return f"{value:.2g}"


def style_log_axis(ax: plt.Axes) -> None:
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(plain_log_tick))
    ax.grid(True, which="major", alpha=0.25)
    ax.grid(True, which="minor", alpha=0.08)


def plot_runtime(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 7.8), sharex=True)
    for row, problem in enumerate(PROBLEM_LABELS):
        for col, method in enumerate(["weno2", "weno3", "ader3"]):
            ax = axes[row, col]
            sub = df[(df["problem"] == problem) & (df["method"] == method)]
            for backend in ["serial", "openmp", "cuda"]:
                b = sub[sub["backend"] == backend].sort_values("nx")
                if b.empty:
                    continue
                ax.plot(
                    b["nx"],
                    b["median_seconds"],
                    marker=MARKERS[backend],
                    color=COLORS[backend],
                    linewidth=1.6,
                    markersize=4.8,
                    label=BACKEND_LABELS[backend],
                )
            ax.set_title(METHOD_LABELS[method], fontsize=11)
            ax.set_xlim(0, 625)
            ax.set_xticks([0, 100, 200, 300, 400, 500, 600])
            ax.set_xlabel("resolution n")
            style_log_axis(ax)
            if col == 0:
                ax.set_ylabel(f"{PROBLEM_LABELS[problem]}\nmedian runtime (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=True, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "hll_discontinuous_full_refresh_runtime.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_speedup(df: pd.DataFrame) -> Path:
    rows = []
    for keys, group in df.groupby(["problem", "method", "nx"]):
        base = group[group["backend"] == "serial"]
        if base.empty:
            continue
        serial_time = float(base["median_seconds"].iloc[0])
        for backend in ["openmp", "cuda"]:
            b = group[group["backend"] == backend]
            if b.empty:
                continue
            rows.append(
                {
                    "problem": keys[0],
                    "method": keys[1],
                    "nx": keys[2],
                    "backend": backend,
                    "speedup": serial_time / float(b["median_seconds"].iloc[0]),
                }
            )
    speed = pd.DataFrame(rows)

    fig, axes = plt.subplots(2, 3, figsize=(14.4, 7.8), sharex=True)
    for row, problem in enumerate(PROBLEM_LABELS):
        for col, method in enumerate(["weno2", "weno3", "ader3"]):
            ax = axes[row, col]
            sub = speed[(speed["problem"] == problem) & (speed["method"] == method)]
            for backend in ["openmp", "cuda"]:
                b = sub[sub["backend"] == backend].sort_values("nx")
                if b.empty:
                    continue
                ax.plot(
                    b["nx"],
                    b["speedup"],
                    marker=MARKERS[backend],
                    color=COLORS[backend],
                    linewidth=1.6,
                    markersize=4.8,
                    label=BACKEND_LABELS[backend],
                )
            ax.axhline(1.0, color="0.45", linestyle="--", linewidth=1.0)
            ax.set_title(METHOD_LABELS[method], fontsize=11)
            ax.set_xlim(0, 625)
            ax.set_xticks([0, 100, 200, 300, 400, 500, 600])
            ax.set_xlabel("resolution n")
            ax.set_ylim(bottom=0.0)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
            ax.grid(True, which="major", alpha=0.25)
            if col == 0:
                ax.set_ylabel(f"{PROBLEM_LABELS[problem]}\nspeedup vs serial")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=True, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "hll_discontinuous_full_refresh_speedup.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> int:
    df = pd.read_csv(CSV)
    for col in ["nx", "ny", "cfl", "median_seconds", "best_seconds"]:
        df[col] = pd.to_numeric(df[col])
    paths = [plot_runtime(df), plot_speedup(df)]
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
