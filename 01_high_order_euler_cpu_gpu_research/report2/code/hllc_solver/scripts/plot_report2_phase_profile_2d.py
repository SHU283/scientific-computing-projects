#!/usr/bin/env python3
"""Plot absolute cost, cost fractions, and phase speedups from phase profiles."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
BACKENDS = ("serial", "openmp", "cuda")
COLORS = {
    "setup/workspace/transfer": "#7f8c8d",
    "CFL reduction": "#d55e00",
    "spatial operator": "#0072b2",
    "update/validity": "#009e73",
    "other": "#e6ab02",
}


def method_label(method: str) -> str:
    return {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}[method]


def load_rows(path: Path) -> dict[tuple[str, str, int], dict[str, float]]:
    grouped: dict[tuple[str, str, int], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            grouped[(row["method"], row["backend"], int(row["nx"]))].append(row)

    medians: dict[tuple[str, str, int], dict[str, float]] = {}
    fields = (
        "total_s",
        "setup_s",
        "cfl_s",
        "workspace_s",
        "spatial_operator_s",
        "state_update_s",
        "admissibility_s",
        "transfer_s",
        "unattributed_s",
    )
    import statistics

    for key, rows in grouped.items():
        medians[key] = {
            field: statistics.median(float(row[field]) for row in rows)
            for field in fields
        }
    return medians


def phase_values(row: dict[str, float]) -> dict[str, float]:
    return {
        "setup/workspace/transfer": row["setup_s"] + row["transfer_s"] + row["workspace_s"],
        "CFL reduction": row["cfl_s"],
        "spatial operator": row["spatial_operator_s"],
        "update/validity": row["state_update_s"] + row["admissibility_s"],
        "other": row["unattributed_s"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    import numpy as np

    windows_arial = Path(r"C:\Windows\Fonts\arial.ttf")
    if windows_arial.exists():
        font_manager.fontManager.addfont(windows_arial)
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        font_family = "Arial"
    except ValueError:
        font_family = "DejaVu Sans"
    style = {"font.family": font_family, "font.size": 15}
    if font_family == "Arial":
        style.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.fallback": "None",
        })
    plt.rcParams.update(style)
    rows = load_rows(args.csv_path)
    resolutions = sorted({key[2] for key in rows})
    if not resolutions:
        raise SystemExit("no phase-profile rows found")
    output_dir = args.output_dir or args.csv_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    phase_names = tuple(COLORS)

    backend_colors = {"serial": "#2b5b84", "openmp": "#d97706", "cuda": "#14877e"}

    # Absolute phase costs are grouped rather than stacked: stacking on a logarithmic
    # axis makes segment heights visually misleading.
    fig, axes = plt.subplots(
        len(resolutions), len(METHODS), figsize=(14.5, 4.2 * len(resolutions)), squeeze=False
    )
    x = np.arange(len(phase_names))
    width = 0.25
    for row_index, resolution in enumerate(resolutions):
        for col_index, method in enumerate(METHODS):
            ax = axes[row_index][col_index]
            for backend_index, backend in enumerate(BACKENDS):
                values = [phase_values(rows[(method, backend, resolution)])[phase] for phase in phase_names]
                visible = [value if value > 0.0 else np.nan for value in values]
                offset = (backend_index - 1) * width
                ax.bar(x + offset, visible, width, color=backend_colors[backend], label=backend)
            ax.set_yscale("log")
            ax.set_xticks(x, ("setup", "CFL", "spatial", "update/validity", "other"), rotation=18)
            ax.set_title(f"{method_label(method)} at {resolution} x {resolution}")
            ax.grid(axis="y", alpha=0.22)
            if col_index == 0:
                ax.set_ylabel("accumulated phase time (s)")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=3, frameon=True)
    fig.suptitle("Absolute HLLC phase costs by backend", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    output = output_dir / "phase_profile_2d_hllc.seconds.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {output}")

    # Fractions answer a different question: which phase limits each complete run?
    fig, axes = plt.subplots(
        len(resolutions), len(METHODS), figsize=(13.5, 4.2 * len(resolutions)), squeeze=False
    )
    for row_index, resolution in enumerate(resolutions):
        for col_index, method in enumerate(METHODS):
            ax = axes[row_index][col_index]
            bottoms = np.zeros(len(BACKENDS))
            for phase in phase_names:
                values = [
                    phase_values(rows[(method, backend, resolution)])[phase]
                    / rows[(method, backend, resolution)]["total_s"]
                    for backend in BACKENDS
                ]
                ax.bar(BACKENDS, values, bottom=bottoms, color=COLORS[phase], label=phase)
                bottoms += np.asarray(values)
            ax.set_ylim(0.0, 1.0)
            ax.set_title(f"{method_label(method)} at {resolution} x {resolution}")
            ax.grid(axis="y", alpha=0.22)
            if col_index == 0:
                ax.set_ylabel("fraction of solver time")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=5, frameon=True)
    fig.suptitle("Where each HLLC run spends its time", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    output = output_dir / "phase_profile_2d_hllc.fraction.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {output}")

    aggregate_phases = (
        "setup/workspace/transfer",
        "CFL reduction",
        "spatial operator",
        "update/validity",
    )
    fig, axes = plt.subplots(len(resolutions), len(METHODS), figsize=(14, 4.0 * len(resolutions)), squeeze=False)
    width = 0.36
    x = np.arange(len(aggregate_phases))
    for row_index, resolution in enumerate(resolutions):
        for col_index, method in enumerate(METHODS):
            ax = axes[row_index][col_index]
            serial_phases = phase_values(rows[(method, "serial", resolution)])
            for offset, backend, color in ((-width / 2, "openmp", "#d55e00"), (width / 2, "cuda", "#009e73")):
                backend_phases = phase_values(rows[(method, backend, resolution)])
                speedups = []
                for phase in aggregate_phases:
                    denominator = backend_phases[phase]
                    speedups.append(serial_phases[phase] / denominator if denominator > 0 else float("nan"))
                ax.bar(x + offset, speedups, width, label=backend, color=color)
            ax.axhline(1.0, color="#555555", linestyle="--", linewidth=1)
            ax.set_yscale("log")
            ax.set_xticks(x, ("setup", "CFL", "spatial", "update/validity"), rotation=18)
            ax.set_title(f"{method_label(method)} at {resolution} x {resolution}")
            ax.grid(axis="y", alpha=0.22)
            if col_index == 0:
                ax.set_ylabel("phase speedup vs serial")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=2, frameon=True)
    fig.suptitle("Where OpenMP and CUDA accelerate the HLLC solver", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    output = output_dir / "phase_profile_2d_hllc.speedup.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {output}")

    # Compact report candidate: compare phase shares across every method/backend
    # without requiring six separate subplot legends.
    fig, axes = plt.subplots(len(resolutions), 1, figsize=(19, 11.6), squeeze=False)
    bar_labels = [
        f"{method_label(method)}\n{backend}"
        for method in METHODS
        for backend in BACKENDS
    ]
    for row_index, resolution in enumerate(resolutions):
        ax = axes[row_index][0]
        bottoms = np.zeros(len(bar_labels))
        for phase in phase_names:
            values = []
            for method in METHODS:
                for backend in BACKENDS:
                    row = rows[(method, backend, resolution)]
                    values.append(phase_values(row)[phase] / row["total_s"])
            bars = ax.bar(
                np.arange(len(bar_labels)), values, bottom=bottoms,
                color=COLORS[phase], label=phase,
            )
            for bar, value, bottom in zip(bars, values, bottoms):
                if value >= 0.08:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bottom + value / 2,
                        f"{100.0 * value:.0f}%",
                        ha="center", va="center", fontsize=14, fontweight="semibold",
                        color="white" if phase == "spatial operator" else "#222222",
                    )
            bottoms += np.asarray(values)
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks(np.arange(len(bar_labels)), bar_labels, fontsize=15)
        ax.tick_params(axis="y", labelsize=15)
        ax.set_ylabel("fraction of profiled solver time", fontsize=17)
        ax.set_title(f"{resolution} x {resolution}", fontsize=19, pad=14)
        ax.grid(axis="y", alpha=0.22)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
        for boundary in (2.5, 5.5):
            ax.axvline(boundary, color="#bbbbbb", linewidth=0.8)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=5,
        frameon=True,
        fontsize=15,
    )
    fig.suptitle("HLLC solver phase breakdown", y=0.995, fontsize=24)
    fig.tight_layout(rect=(0.025, 0.02, 0.99, 0.89), h_pad=2.4)
    output = output_dir / "phase_profile_2d_hllc.breakdown.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {output}")

    # Focused large-grid figure for the performance chapter.  Phase speedups are
    # relative to serial for the same method at n=600.
    large_resolution = max(resolutions)
    fig, axes = plt.subplots(
        1,
        len(METHODS),
        figsize=(17, 9.6),
        squeeze=False,
    )
    comparable_phases = (
        "CFL reduction",
        "spatial operator",
        "update/validity",
    )
    phase_labels = ("CFL\nreduction", "spatial\noperator", "update /\nvalidity")
    x = np.arange(len(comparable_phases))
    width = 0.34
    for col_index, method in enumerate(METHODS):
        ax = axes[0][col_index]
        serial_phases = phase_values(rows[(method, "serial", large_resolution)])
        panel_max = 1.0
        for offset, backend, color in (
            (-width / 2, "openmp", "#d97706"),
            (width / 2, "cuda", "#14877e"),
        ):
            backend_phases = phase_values(rows[(method, backend, large_resolution)])
            speedups = [
                serial_phases[phase] / backend_phases[phase]
                if backend_phases[phase] > 0.0 else float("nan")
                for phase in comparable_phases
            ]
            finite_speedups = [value for value in speedups if np.isfinite(value)]
            if finite_speedups:
                panel_max = max(panel_max, max(finite_speedups))
            bars = ax.bar(
                x + offset,
                speedups,
                width=width,
                label=backend,
                color=color,
            )
            for bar, value in zip(bars, speedups):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value * 1.06,
                    f"{value:.1f}x" if value < 100.0 else f"{value:.0f}x",
                    ha="center",
                    va="bottom",
                    fontsize=13,
                )
        ax.axhline(1.0, color="#555555", linestyle="--", linewidth=1)
        ax.set_yscale("log")
        ax.set_ylim(0.8, panel_max * 1.55)
        ax.set_xticks(x, phase_labels)
        ax.set_title(method_label(method), fontsize=18, pad=11)
        ax.set_xlabel("solver phase", fontsize=16, labelpad=8)
        ax.tick_params(axis="both", labelsize=15)
        ax.grid(axis="y", alpha=0.22)
        if col_index == 0:
            ax.set_ylabel("phase speedup relative to serial", fontsize=16)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=2,
        frameon=True,
        fontsize=15,
    )
    fig.suptitle(
        f"Comparable phase speedups at {large_resolution} x {large_resolution}",
        y=0.995, fontsize=21,
    )
    fig.tight_layout(rect=(0.015, 0.02, 0.99, 0.88), w_pad=2.8)
    output = output_dir / "phase_profile_2d_hllc.n600_phase_speedup.png"
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
