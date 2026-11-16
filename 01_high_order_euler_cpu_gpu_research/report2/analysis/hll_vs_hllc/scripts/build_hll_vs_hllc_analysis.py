from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd


OUT = Path(__file__).resolve().parents[1]
PROCESSED = OUT / "processed_results"
FIGURES = OUT / "figures" / "hll_vs_hllc"
SOURCE = OUT / "source_data"

HLL_SMOOTH = SOURCE / "hll" / "smooth.csv"
HLLC_SMOOTH = SOURCE / "hllc" / "smooth.csv"

HLL_DISCONT = SOURCE / "hll" / "discontinuous.csv"
HLLC_DISCONT = SOURCE / "hllc" / "discontinuous.csv"

HLL_GPU = SOURCE / "hll" / "gpu_limit.csv"
HLLC_GPU = SOURCE / "hllc" / "gpu_limit.csv"

HLL_OMP = SOURCE / "hll" / "openmp_thread_scaling.csv"
HLLC_OMP = SOURCE / "hllc" / "openmp_thread_scaling.csv"


COLORS = {"HLL": "#1f5a8a", "HLLC": "#0c8a7a"}
METHODS = ["weno2", "weno3", "ader3"]
METHOD_LABELS = {"weno2": "WENO-2", "weno3": "WENO-3", "ader3": "ADER-3"}
BACKEND_MARKERS = {"serial": "o", "openmp": "s", "cuda": "^"}


def ensure_dirs() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path, flux: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    for col in ["problem", "method", "backend"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip().str.strip('"').str.lower()
    for col in df.columns:
        if col not in {"problem", "method", "backend", "samples_seconds"}:
            converted = pd.to_numeric(df[col], errors="coerce")
            if converted.notna().any():
                df[col] = converted
    df.insert(0, "flux", flux)
    return df


def savefig(name: str) -> Path:
    path = FIGURES / name
    fig = plt.gcf()
    for ax in fig.axes:
        if ax.get_xscale() == "log" and not getattr(ax, "_preserve_x_formatter", False):
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(plain_sci_tick))
        if ax.get_yscale() == "log" and not getattr(ax, "_preserve_y_formatter", False):
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(plain_sci_tick))
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()
    return path


def plain_sci_tick(value: float, _pos: int | None = None) -> str:
    if value <= 0 or not math.isfinite(value):
        return ""
    exponent = round(math.log10(value))
    if math.isclose(value, 10**exponent, rel_tol=1e-9, abs_tol=0.0):
        if exponent == 0:
            return "1"
        return f"1e{exponent}"
    return f"{value:g}"


def compact_runtime_tick(value: float, _pos: int | None = None) -> str:
    if value <= 0 or not math.isfinite(value):
        return ""
    if value < 1:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    if value < 10:
        return f"{value:.0f}" if float(value).is_integer() else f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.0f}"


def style_log_axis(
    ax: plt.Axes,
    axis: str,
    *,
    formatter=plain_sci_tick,
    subs: tuple[float, ...] = (1.0, 2.0, 5.0),
) -> None:
    """Use readable 1/2/5 log ticks instead of only decade ticks."""
    target = ax.xaxis if axis == "x" else ax.yaxis
    target.set_major_locator(ticker.LogLocator(base=10, subs=subs))
    target.set_major_formatter(ticker.FuncFormatter(formatter))
    target.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(1.0, 10.0)))
    target.set_minor_formatter(ticker.NullFormatter())
    if axis == "x":
        ax._preserve_x_formatter = True
    else:
        ax._preserve_y_formatter = True


def plot_style() -> None:
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in ["Arial", "Calibri", "Segoe UI", "Liberation Sans"]:
        if candidate in available_fonts:
            font_family = candidate
            break
    else:
        font_family = sorted(available_fonts)[0]
    plt.rcParams.update(
        {
            "font.family": font_family,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.titlesize": 18,
            "axes.labelsize": 16,
            "legend.fontsize": 15,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "figure.dpi": 120,
            "axes.formatter.use_mathtext": False,
            "text.usetex": False,
        }
    )


def exact_match(left: pd.DataFrame, right: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    a = left.copy()
    b = right.copy()
    suffixes = ("_hll", "_hllc")
    return a.merge(b, on=keys, how="inner", suffixes=suffixes)


def practical_match(left: pd.DataFrame, right: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Match runs that share the same physical case, method, backend, grid, and final time.

    CFL is intentionally not part of these keys.  This is used only for the practical
    HLL-vs-HLLC discontinuous comparison where ADER-3 uses route-specific stable CFL
    values.
    """
    a = left.copy()
    b = right.copy()
    suffixes = ("_hll", "_hllc")
    return a.merge(b, on=keys, how="inner", suffixes=suffixes)


def add_speedups(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    key_cols = [c for c in ["flux", "problem", "method", "nx", "ny", "cfl", "final_time"] if c in out.columns]
    serial = out[out["backend"] == "serial"][key_cols + ["best_seconds"]].rename(columns={"best_seconds": "serial_seconds"})
    out = out.merge(serial, on=key_cols, how="left")
    out["speedup_vs_serial"] = out["serial_seconds"] / out["best_seconds"]
    return out


def smooth_analysis(hll: pd.DataFrame, hllc: pd.DataFrame) -> dict[str, pd.DataFrame]:
    keys = ["problem", "method", "backend", "nx", "ny", "cfl", "final_time"]
    matched = exact_match(hll, hllc, keys)
    matched.to_csv(PROCESSED / "smooth_matched_hll_vs_hllc.csv", index=False)

    n600 = matched[matched["nx"] == 600].copy()
    n600.to_csv(PROCESSED / "smooth_n600_matched_hll_vs_hllc.csv", index=False)

    # Errors are independent of backend for these smooth runs; use serial as the clean representative.
    serial_errors = n600[n600["backend"] == "serial"].copy()
    serial_errors.to_csv(PROCESSED / "smooth_n600_error_norms_hll_vs_hllc.csv", index=False)

    rel = matched[matched["backend"] == "serial"].copy()
    for metric in ["l1_rho", "l2_rho", "linf_rho"]:
        rel[f"{metric}_relative_difference"] = (rel[f"{metric}_hllc"] - rel[f"{metric}_hll"]) / rel[f"{metric}_hll"]
    rel.to_csv(PROCESSED / "smooth_error_relative_difference.csv", index=False)

    return {"matched": matched, "n600": n600, "serial_errors": serial_errors, "relative": rel}


def plot_smooth_candidates(data: dict[str, pd.DataFrame]) -> list[Path]:
    paths: list[Path] = []
    n600 = data["n600"]
    serial_errors = data["serial_errors"]
    rel = data["relative"]
    matched = data["matched"]

    metrics = [("l1_rho", "L1 density error"), ("l2_rho", "L2 density error"), ("linf_rho", "Linf density error")]
    x = np.arange(len(METHODS))
    width = 0.36
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True)
    for ax, (metric, label) in zip(axes, metrics):
        hll_vals = [serial_errors.loc[serial_errors["method"] == m, f"{metric}_hll"].iloc[0] for m in METHODS]
        hllc_vals = [serial_errors.loc[serial_errors["method"] == m, f"{metric}_hllc"].iloc[0] for m in METHODS]
        ax.bar(x - width / 2, hll_vals, width, label="HLL", color=COLORS["HLL"])
        ax.bar(x + width / 2, hllc_vals, width, label="HLLC", color=COLORS["HLLC"])
        ax.set_title(label)
        ax.set_yscale("log")
        style_log_axis(ax, "y")
        ax.set_xticks(x, [METHOD_LABELS[m] for m in METHODS])
        ax.set_ylabel("error")
    axes[0].legend(loc="best")
    fig.suptitle("Smooth n=600: HLL vs HLLC density error norms")
    paths.append(savefig("smooth_n600_error_norms.png"))

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.6), sharey=True)
    for ax, method in zip(axes, METHODS):
        subset = n600[n600["method"] == method].copy()
        for flux, color in COLORS.items():
            suffix = "_hll" if flux == "HLL" else "_hllc"
            xs = []
            ys = []
            for backend in ["serial", "openmp", "cuda"]:
                row = subset[subset["backend"] == backend].iloc[0]
                xval = row[f"best_seconds{suffix}"]
                yval = row[f"l1_rho{suffix}"]
                xs.append(xval)
                ys.append(yval)
                ax.scatter(
                    xval,
                    yval,
                    color=color,
                    marker=BACKEND_MARKERS.get(backend, "o"),
                    s=76,
                    zorder=3,
                )
            ax.plot(xs, ys, color=color, linewidth=1.5, alpha=0.55)
        ax.set_xscale("log")
        ax.set_yscale("log")
        style_log_axis(ax, "x", formatter=compact_runtime_tick)
        style_log_axis(ax, "y")
        ax.set_xlabel("best runtime (s)", fontsize=16)
        ax.set_title(METHOD_LABELS[method], fontsize=18, pad=11)
        ax.tick_params(axis="both", which="major", labelsize=15)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
    axes[0].set_ylabel("L1 density error", fontsize=16)
    for backend, marker in BACKEND_MARKERS.items():
        axes[0].scatter([], [], color="gray", marker=marker, label=backend)
    for flux, color in COLORS.items():
        axes[0].plot([], [], color=color, marker="o", label=flux)
    axes[0].legend(ncol=1, loc="best", fontsize=15)
    fig.suptitle("Smooth n=600 accuracy-runtime trade-off", fontsize=21, y=0.99)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.84, wspace=0.20)
    paths.append(savefig("smooth_n600_runtime_vs_l1.png"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    serial = matched[matched["backend"] == "serial"]
    for ax, method in zip(axes, METHODS):
        sub = serial[serial["method"] == method].sort_values("nx")
        xpos = np.arange(len(sub))
        ax.plot(xpos, sub["l1_rho_hll"], marker="o", color=COLORS["HLL"], label="HLL")
        ax.plot(xpos, sub["l1_rho_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
        ax.set_yscale("log")
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlabel("resolution n")
        ax.set_xticks(xpos, [str(int(n)) for n in sub["nx"]])
    axes[0].set_ylabel("L1 density error")
    axes[0].legend()
    fig.suptitle("Smooth convergence comparison")
    paths.append(savefig("smooth_convergence_hll_vs_hllc.png"))

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    for method in METHODS:
        sub = rel[rel["method"] == method].sort_values("nx")
        ax.plot(sub["nx"], 100.0 * sub["l1_rho_relative_difference"], marker="o", label=METHOD_LABELS[method])
    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("resolution n")
    ax.set_ylabel("(HLLC L1 - HLL L1) / HLL L1 (%)")
    ax.set_title("Relative L1-density error difference")
    ax.legend()
    paths.append(savefig("smooth_relative_l1_difference.png"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)
    for ax, method in zip(axes, METHODS):
        sub = matched[(matched["backend"] == "cuda") & (matched["method"] == method)].sort_values("nx")
        ax.plot(sub["nx"], sub["best_seconds_hll"], marker="o", color=COLORS["HLL"], label="HLL")
        ax.plot(sub["nx"], sub["best_seconds_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
        ax.set_yscale("log")
        style_log_axis(ax, "y", formatter=compact_runtime_tick)
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlabel("resolution n")
        ax.set_ylabel("CUDA best runtime (s)")
    axes[0].legend()
    fig.suptitle("Smooth CUDA runtime comparison")
    paths.append(savefig("smooth_cuda_runtime_hll_vs_hllc.png"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, method in zip(axes, METHODS):
        sub = matched[(matched["backend"] == "cuda") & (matched["method"] == method)].sort_values("nx")
        ax.plot(
            sub["nx"],
            sub["l1_rho_hllc"],
            marker="s",
            color=COLORS["HLLC"],
            linewidth=1.6,
            label="HLLC",
            zorder=2,
        )
        ax.plot(
            sub["nx"],
            sub["l1_rho_hll"],
            marker="o",
            markerfacecolor="white",
            markeredgewidth=1.2,
            color=COLORS["HLL"],
            linestyle="--",
            linewidth=1.4,
            label="HLL",
            zorder=3,
        )
        ax.set_yscale("log")
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlabel("resolution n")
    axes[0].set_ylabel("L1 density error")
    axes[0].legend()
    fig.suptitle("Smooth CUDA density-error comparison")
    paths.append(savefig("smooth_cuda_l1_hll_vs_hllc.png"))

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.6), sharey=True)
    for ax, method in zip(axes, METHODS):
        sub = matched[(matched["backend"] == "cuda") & (matched["method"] == method)].sort_values("nx")
        ax.plot(
            sub["best_seconds_hll"],
            sub["l1_rho_hll"],
            marker="o",
            markersize=7,
            linewidth=1.8,
            color=COLORS["HLL"],
            label="HLL",
        )
        ax.plot(
            sub["best_seconds_hllc"],
            sub["l1_rho_hllc"],
            marker="s",
            markersize=7,
            linewidth=1.8,
            color=COLORS["HLLC"],
            label="HLLC",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        style_log_axis(ax, "x", formatter=compact_runtime_tick)
        style_log_axis(ax, "y", subs=(1.0,))
        ax.set_title(METHOD_LABELS[method], fontsize=18, pad=11)
        ax.set_xlabel("CUDA best runtime (s)", fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=15)
        for spine in ax.spines.values():
            spine.set_linewidth(1.1)
    axes[0].set_ylabel("L1 density error", fontsize=16)
    axes[0].legend(fontsize=15, loc="best")
    fig.suptitle("Smooth CUDA accuracy-runtime comparison", fontsize=21, y=0.99)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.84, wspace=0.20)
    paths.append(savefig("smooth_cuda_runtime_vs_l1_hll_vs_hllc.png"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    backends = ["serial", "openmp", "cuda"]
    xpos = np.arange(len(backends))
    for ax, method in zip(axes, METHODS):
        sub = n600[n600["method"] == method].set_index("backend")
        ax.bar(xpos - width / 2, [sub.loc[b, "best_seconds_hll"] for b in backends], width, label="HLL", color=COLORS["HLL"])
        ax.bar(xpos + width / 2, [sub.loc[b, "best_seconds_hllc"] for b in backends], width, label="HLLC", color=COLORS["HLLC"])
        ax.set_xticks(xpos, backends)
        ax.set_yscale("log")
        ax.set_title(METHOD_LABELS[method])
        ax.set_ylabel("best runtime (s)")
    axes[0].legend()
    fig.suptitle("Smooth n=600 runtime comparison")
    paths.append(savefig("smooth_n600_runtime_hll_vs_hllc.png"))

    return paths


def runtime_analysis(hll: pd.DataFrame, hllc: pd.DataFrame, name: str, keys: list[str]) -> pd.DataFrame:
    matched = exact_match(hll, hllc, keys)
    matched.to_csv(PROCESSED / f"{name}_matched_hll_vs_hllc.csv", index=False)
    return matched


def practical_runtime_analysis(hll: pd.DataFrame, hllc: pd.DataFrame, name: str, keys: list[str]) -> pd.DataFrame:
    matched = practical_match(hll, hllc, keys)
    matched.to_csv(PROCESSED / f"{name}_practical_hll_vs_hllc.csv", index=False)
    return matched


def plot_gpu_limit(gpu_match: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if gpu_match.empty:
        return paths
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.4), sharey=False)
    for ax, method in zip(axes, METHODS):
        sub = gpu_match[gpu_match["method"] == method].sort_values("nx")
        ax.plot(sub["nx"], sub["best_seconds_hll"], marker="o", color=COLORS["HLL"], label="HLL")
        ax.plot(sub["nx"], sub["best_seconds_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
        ax.set_yscale("log")
        style_log_axis(ax, "y", formatter=compact_runtime_tick)
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlim(450, 5150)
        ax.set_xticks([600, 1800, 3000, 4200, 5000])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.91))
    fig.suptitle("CUDA-only large-grid runtime comparison", y=0.985, fontsize=21)
    fig.supxlabel("resolution n", y=0.04, fontsize=18)
    fig.supylabel("CUDA runtime (s)", x=0.018, fontsize=18)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.76, wspace=0.20)
    paths.append(savefig("gpu_limit_cuda_runtime_hll_vs_hllc.png"))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, method in zip(axes, METHODS):
        sub = gpu_match[gpu_match["method"] == method].sort_values("nx")
        ax.plot(sub["nx"], sub["l1_rho_hll"], marker="o", color=COLORS["HLL"], label="HLL")
        ax.plot(sub["nx"], sub["l1_rho_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlabel("resolution n")
    axes[0].set_ylabel("L1 density error")
    axes[0].legend()
    fig.suptitle("CUDA-only large-grid error comparison")
    paths.append(savefig("gpu_limit_l1_hll_vs_hllc.png"))
    return paths


def plot_discontinuous(dis_match: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if dis_match.empty:
        return paths
    cuda = dis_match[dis_match["backend"] == "cuda"].copy()
    if not cuda.empty:
        problems = ["quadrant_riemann_2d", "shock_bubble_2d"]
        methods = [method for method in METHODS if (cuda["method"] == method).any()]
        problem_labels = {
            "quadrant_riemann_2d": "Quadrant Riemann",
            "shock_bubble_2d": "Shock-bubble",
        }
        fig, axes = plt.subplots(
            len(problems),
            len(methods),
            figsize=(12.8, 8.0),
            sharex=True,
            sharey=True,
        )
        axes = np.atleast_2d(axes)
        for i, problem in enumerate(problems):
            for j, method in enumerate(methods):
                ax = axes[i, j]
                sub = cuda[(cuda["problem"] == problem) & (cuda["method"] == method)].sort_values("nx")
                if sub.empty:
                    ax.axis("off")
                    continue
                ax.plot(sub["nx"], sub["median_seconds_hll"], marker="o", color=COLORS["HLL"], label="HLL")
                ax.plot(sub["nx"], sub["median_seconds_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
                ax.set_yscale("log")
                style_log_axis(ax, "y", formatter=compact_runtime_tick)
                ax.set_xlim(25, 625)
                ax.set_xticks([50, 100, 200, 400, 600])
                ax.set_title(f"{problem_labels[problem]} / {METHOD_LABELS[method]}")

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.92))
        fig.suptitle("Discontinuous matched CUDA runtime comparison", y=0.985, fontsize=20)
        fig.supxlabel("resolution n", y=0.025, fontsize=17)
        fig.supylabel("CUDA median runtime (s)", x=0.02, fontsize=17)
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.11, top=0.80, wspace=0.10, hspace=0.28)
        paths.append(savefig("discontinuous_cuda_runtime_matched.png"))
    return paths


def plot_discontinuous_practical(dis_practical: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if dis_practical.empty:
        return paths
    cuda = dis_practical[dis_practical["backend"] == "cuda"].copy()
    if cuda.empty:
        return paths

    problems = ["quadrant_riemann_2d", "shock_bubble_2d"]
    problem_labels = {
        "quadrant_riemann_2d": "quadrant riemann",
        "shock_bubble_2d": "shock bubble",
    }
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.2), sharex=True)
    for i, problem in enumerate(problems):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = cuda[(cuda["problem"] == problem) & (cuda["method"] == method)].sort_values("nx")
            if sub.empty:
                ax.axis("off")
                continue
            ax.plot(sub["nx"], sub["median_seconds_hll"], marker="o", color=COLORS["HLL"], label="HLL")
            ax.plot(sub["nx"], sub["median_seconds_hllc"], marker="s", color=COLORS["HLLC"], label="HLLC")
            ax.set_yscale("log")
            style_log_axis(ax, "y", formatter=compact_runtime_tick)
            ax.set_xlim(0, 625)
            ax.set_xticks([0, 100, 200, 300, 400, 500, 600])
            ax.set_title(f"{problem_labels[problem]} / {METHOD_LABELS[method]}")
            ax.set_xlabel("resolution n")
            if j == 0:
                ax.set_ylabel("CUDA median runtime (s)")
            if method == "ader3":
                cfl_pairs = sorted(
                    {
                        (float(row["cfl_hll"]), float(row["cfl_hllc"]))
                        for _, row in sub.iterrows()
                        if "cfl_hll" in sub.columns and "cfl_hllc" in sub.columns
                    }
                )
                if cfl_pairs:
                    ax.text(
                        0.03,
                        0.95,
                        "stable CFLs:\n"
                        + "\n".join(f"HLL {a:g}, HLLC {b:g}" for a, b in cfl_pairs),
                        transform=ax.transAxes,
                        va="top",
                        ha="left",
                        fontsize=7,
                        bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.75},
                    )
    axes[0, 0].legend()
    fig.suptitle("Discontinuous CUDA runtime comparison, practical stable-CFL view")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    paths.append(savefig("discontinuous_cuda_runtime_practical_hll_vs_hllc.png"))
    return paths


def plot_openmp(omp_match: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if omp_match.empty:
        return paths
    rows = []
    keys = ["problem", "backend", "method", "nx", "ny", "cfl", "final_time", "threads"]
    for _, row in omp_match.iterrows():
        base = {k: row[k] for k in keys}
        for flux, suffix in [("HLL", "_hll"), ("HLLC", "_hllc")]:
            rec = base.copy()
            rec["flux"] = flux
            rec["best_seconds"] = row[f"best_seconds{suffix}"]
            rows.append(rec)
    both = pd.DataFrame(rows)
    base_keys = ["flux", "problem", "method", "nx", "ny", "cfl", "final_time"]
    thread_one = (
        both[both["threads"] == 1][base_keys + ["best_seconds"]]
        .rename(columns={"best_seconds": "thread_one_seconds"})
    )
    both = both.merge(thread_one, on=base_keys, how="left")
    both["speedup_vs_thread_one"] = both["thread_one_seconds"] / both["best_seconds"]
    both.to_csv(PROCESSED / "openmp_thread_matched_speedups_hll_vs_hllc.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.4), sharex=True, sharey=True)
    for ax, method in zip(axes, METHODS):
        sub = both[both["method"] == method].sort_values("threads")
        for flux, color in COLORS.items():
            s = sub[sub["flux"] == flux]
            if not s.empty:
                ax.plot(
                    s["threads"],
                    s["speedup_vs_thread_one"],
                    marker="o" if flux == "HLL" else "s",
                    linewidth=2.0,
                    color=color,
                    label=flux,
                )
        ax.plot([1, 2, 4, 8, 16], [1, 2, 4, 8, 16], color="gray", linestyle="--", linewidth=1, label="ideal")
        ax.set_title(METHOD_LABELS[method])
        ax.set_xticks([1, 2, 4, 8, 16])
        ax.set_xlim(0.5, 16.5)
        ax.set_ylim(0.5, 16.8)
        ax.set_yticks([1, 2, 4, 8, 12, 16])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.91))
    fig.suptitle("OpenMP thread scaling: HLL vs HLLC", y=0.985, fontsize=21)
    fig.supxlabel("OpenMP threads", y=0.04, fontsize=18)
    fig.supylabel("speedup relative to 1 thread", x=0.018, fontsize=18)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.76, wspace=0.10)
    paths.append(savefig("openmp_thread_speedup_hll_vs_hllc.png"))
    return paths


def write_summary(
    smooth: dict[str, pd.DataFrame],
    dis: pd.DataFrame,
    dis_practical: pd.DataFrame,
    gpu: pd.DataFrame,
    omp: pd.DataFrame,
    figures: Iterable[Path],
) -> None:
    figure_lines = [f"- `{p.relative_to(OUT).as_posix()}`" for p in figures]
    smooth_match = smooth["matched"]
    n600 = smooth["n600"]
    rel = smooth["relative"]

    dis_exact = len(dis)
    dis_methods = ", ".join(sorted(dis["method"].unique())) if dis_exact else "none"
    dis_practical_count = len(dis_practical)
    dis_practical_methods = ", ".join(sorted(dis_practical["method"].unique())) if dis_practical_count else "none"
    gpu_ns = ", ".join(str(int(x)) for x in sorted(gpu["nx"].unique())) if not gpu.empty else "none"
    omp_methods = ", ".join(sorted(omp["method"].unique())) if not omp.empty else "none"

    rel600 = rel[rel["nx"] == 600][["method", "l1_rho_relative_difference"]].copy()
    rel600["percent"] = 100.0 * rel600["l1_rho_relative_difference"]
    rel_text = "\n".join(
        f"- {METHOD_LABELS[row.method]}: HLLC L1 is {row.percent:+.1f}% relative to HLL at n=600."
        for row in rel600.itertuples()
    )

    text = f"""# HLL vs HLLC Analysis Summary

Generated by `analysis_outputs/hll_vs_hllc/scripts/build_hll_vs_hllc_analysis.py`.

## Matched Experiments Found

- Smooth 2D Gaussian performance: {len(smooth_match)} exact matched rows. Keys matched on problem, method, backend, grid, CFL, and final time.
- Smooth 2D Gaussian n=600 subset: {len(n600)} matched rows.
- Discontinuous runtime: {dis_exact} exact matched rows. Matched methods: {dis_methods}. Rows with different CFL choices are excluded from strict comparison.
- Discontinuous practical stable-CFL runtime: {dis_practical_count} matched rows. Matched methods: {dis_practical_methods}. Here CFL is allowed to differ when each route uses its selected stable value, so ADER-3 is included with a caveat.
- CUDA-only GPU-limit smooth sweep: {len(gpu)} exact matched rows. Matched resolutions: {gpu_ns}.
- OpenMP thread scaling: {len(omp)} exact matched rows. Matched methods: {omp_methods}.

## Metrics Available

- Smooth Gaussian: `best_seconds`, `mean_seconds`, `median_seconds`, `l1_rho`, `l2_rho`, and `linf_rho`.
- GPU-limit smooth sweep: CUDA runtime and density error metrics.
- Discontinuous cases: runtime and speedup information only; no exact density-error norms are recorded in the timing CSVs.
- OpenMP thread scaling: runtime and derived speedup; no accuracy metrics.

## Smooth n=600 Accuracy Notes

{rel_text}

These values should be interpreted from the data, not as an assumption that HLLC is always more accurate.

## Candidate Plots Generated

{chr(10).join(figure_lines)}

## Best 2-4 Plots To Include

1. `figures/hll_vs_hllc/smooth_n600_runtime_vs_l1.png`: most directly supports the Section 3 accuracy-efficiency comparison.
2. `figures/hll_vs_hllc/smooth_n600_error_norms.png`: cleanest accuracy-only comparison at the matched 600-grid point.
3. `figures/hll_vs_hllc/smooth_cuda_runtime_hll_vs_hllc.png`: shows how the CUDA cost changes across grid size for each method.
4. `figures/hll_vs_hllc/gpu_limit_cuda_runtime_hll_vs_hllc.png`: useful if Section 3 discusses large-grid GPU behaviour.

## Technically Correct But Probably Secondary

- `smooth_relative_l1_difference.png`: useful for internal interpretation, but the signed percentage can distract unless the text explains it carefully.
- `smooth_convergence_hll_vs_hllc.png`: useful if convergence is discussed, but the main report can keep HLLC as the validation route and use HLL only in Section 3.
- `discontinuous_cuda_runtime_matched.png`: helpful as strict supporting evidence, but it excludes ADER-3 when CFL choices differ.
- `discontinuous_cuda_runtime_practical_hll_vs_hllc.png`: useful for the practical stable-CFL comparison because it includes ADER-3; label it clearly as not a strict same-CFL comparison.
- `openmp_thread_speedup_hll_vs_hllc.png`: useful if thread scaling is included, but less central than accuracy-efficiency.

## Missing Or Limited Metrics

- No discontinuous accuracy norms are available in the runtime tables, so discontinuous HLL-vs-HLLC accuracy comparisons would require new post-processing against a reference or a field-to-field diagnostic.
- Smooth accuracy metrics currently cover density only (`rho`). Momentum, velocity, pressure, and energy error comparisons would require recomputing diagnostics from field outputs.
- Backend-consistency tables measure backend agreement within one flux route, not HLL-vs-HLLC physical accuracy.
- ADER-3 discontinuous comparisons are excluded from exact same-CFL matching when CFL differs between routes.  They can still be shown in the practical stable-CFL view if the report text states that each route uses its most stable selected CFL.

## Report Integration Rule

HLLC should remain the main Report 2 route for validation, field figures, backend consistency, and conclusions. HLL should appear only as a Section 3 comparator in the accuracy-efficiency and GPU-scaling discussion, and only where the matched-experiment conditions are explicit.
"""
    (OUT / "summary.md").write_text(text, encoding="utf-8")


def main() -> int:
    ensure_dirs()
    plot_style()

    hll_smooth = read_csv(HLL_SMOOTH, "HLL")
    hllc_smooth = read_csv(HLLC_SMOOTH, "HLLC")
    smooth = smooth_analysis(hll_smooth, hllc_smooth)
    figures = plot_smooth_candidates(smooth)

    hll_dis = read_csv(HLL_DISCONT, "HLL")
    hllc_dis = read_csv(HLLC_DISCONT, "HLLC")
    dis_keys = ["problem", "method", "backend", "nx", "ny", "cfl", "final_time"]
    dis_match = runtime_analysis(hll_dis, hllc_dis, "discontinuous_runtime", dis_keys)
    figures.extend(plot_discontinuous(dis_match))
    dis_practical_keys = ["problem", "method", "backend", "nx", "ny", "final_time"]
    dis_practical = practical_runtime_analysis(
        hll_dis,
        hllc_dis,
        "discontinuous_runtime",
        dis_practical_keys,
    )
    figures.extend(plot_discontinuous_practical(dis_practical))

    hll_gpu = read_csv(HLL_GPU, "HLL")
    hllc_gpu = read_csv(HLLC_GPU, "HLLC")
    gpu_keys = ["problem", "method", "backend", "nx", "ny", "cfl", "final_time"]
    gpu_match = runtime_analysis(hll_gpu, hllc_gpu, "gpu_limit", gpu_keys)
    figures.extend(plot_gpu_limit(gpu_match))

    hll_omp = read_csv(HLL_OMP, "HLL")
    hllc_omp = read_csv(HLLC_OMP, "HLLC")
    omp_keys = ["problem", "backend", "method", "nx", "ny", "cfl", "final_time", "threads"]
    omp_match = runtime_analysis(hll_omp, hllc_omp, "openmp_thread", omp_keys)
    figures.extend(plot_openmp(omp_match))

    write_summary(smooth, dis_match, dis_practical, gpu_match, omp_match, figures)

    print(f"Wrote processed tables to: {PROCESSED}")
    print(f"Wrote figures to: {FIGURES}")
    print(f"Wrote summary to: {OUT / 'summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
