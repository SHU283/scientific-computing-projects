#!/usr/bin/env python3
"""Run isolated Report 2 convergence checks for the 1D high-order methods.

Outputs are written below results/report2 so that Report 1 data remains
untouched. Errors are measured against exact finite-volume cell averages.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
PROBLEMS = ("smooth_advection", "gaussian_advection")
DEFAULT_RESOLUTIONS = (50, 100, 200, 400, 800)
METHOD_CFL = {
    "weno2": 0.4,
    "weno3": 0.4,
    "ader3": 0.25,
}
EXPECTED_FINAL_ORDER = {
    "weno2": (1.8, 2.2),
    "weno3": (2.8, 3.2),
    "ader3": (2.8, 3.2),
}


def smooth_average(
    x_left: float,
    x_right: float,
    x_min: float,
    x_max: float,
    phase_shift: float,
) -> float:
    length = x_max - x_min
    wave_number = 2.0 * math.pi / length
    integral = (
        -math.cos(wave_number * (x_right - x_min - phase_shift))
        + math.cos(wave_number * (x_left - x_min - phase_shift))
    ) / wave_number
    return 1.0 + 0.2 * integral / (x_right - x_left)


def gaussian_average(
    x_left: float,
    x_right: float,
    center: float,
    sigma: float,
) -> float:
    sqrt_half = 1.0 / math.sqrt(2.0)
    integral = sigma * math.sqrt(math.pi / 2.0) * (
        math.erf((x_right - center) * sqrt_half / sigma)
        - math.erf((x_left - center) * sqrt_half / sigma)
    )
    return 1.0 + 0.2 * integral / (x_right - x_left)


def periodic_gaussian_average(
    x_left: float,
    x_right: float,
    center: float,
    sigma: float,
    x_min: float,
    x_max: float,
) -> float:
    length = x_max - x_min
    return 1.0 + sum(
        gaussian_average(x_left, x_right, center + shift * length, sigma) - 1.0
        for shift in (-1, 0, 1)
    )


def exact_density(
    problem: str,
    x: float,
    dx: float,
    final_time: float,
    velocity: float,
    x_min: float,
    x_max: float,
) -> float:
    x_left = x - 0.5 * dx
    x_right = x + 0.5 * dx
    length = x_max - x_min

    if problem == "smooth_advection":
        phase_shift = 0.173 * length + velocity * final_time
        return smooth_average(x_left, x_right, x_min, x_max, phase_shift)

    center = x_min + 0.317 * length + velocity * final_time
    center = x_min + ((center - x_min) % length)
    sigma = 0.08 * length
    return periodic_gaussian_average(
        x_left, x_right, center, sigma, x_min, x_max
    )


def load_density(path: Path) -> tuple[list[float], list[float]]:
    with path.open("r", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return (
        [float(row["x"]) for row in rows],
        [float(row["rho"]) for row in rows],
    )


def error_norms(
    numerical: list[float], exact: list[float]
) -> tuple[float, float, float]:
    errors = [abs(value - reference) for value, reference in zip(numerical, exact)]
    l1 = sum(errors) / len(errors)
    l2 = math.sqrt(sum(error * error for error in errors) / len(errors))
    linf = max(errors)
    return l1, l2, linf


def observed_order(coarse_error: float | None, fine_error: float) -> float | None:
    if coarse_error is None or coarse_error <= 0.0 or fine_error <= 0.0:
        return None
    return math.log(coarse_error / fine_error, 2.0)


def check_asymptotic_orders(
    rows: list[dict[str, str | int | float]],
) -> None:
    failures: list[str] = []
    for problem in PROBLEMS:
        for method in METHODS:
            method_rows = [
                row
                for row in rows
                if row["problem"] == problem and row["method"] == method
            ]
            final_row = method_rows[-1]
            lower, upper = EXPECTED_FINAL_ORDER[method]
            for metric in ("order_l1", "order_l2", "order_linf"):
                order = float(final_row[metric])
                if not lower <= order <= upper:
                    failures.append(
                        f"{problem}/{method}/{metric}={order:.3f}, "
                        f"expected [{lower:.1f}, {upper:.1f}]"
                    )

    if failures:
        raise RuntimeError(
            "Report 2 convergence checks failed:\n  " + "\n  ".join(failures)
        )
    print("Report 2 convergence checks: PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2")
    parser.add_argument("--output-dir", default="results/report2/validation_1d")
    parser.add_argument(
        "--resolutions", nargs="+", type=int, default=DEFAULT_RESOLUTIONS
    )
    parser.add_argument("--final-time", type=float, default=0.3)
    parser.add_argument(
        "--backend", choices=("serial", "openmp"), default="serial"
    )
    parser.add_argument(
        "--cfl",
        type=float,
        default=None,
        help="Optional common-CFL override for diagnostics only.",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    executable = root / args.build_dir / "final_project"
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    x_min = 0.0
    x_max = 1.0
    velocity = 0.73
    summary_rows: list[dict[str, str | int | float]] = []

    for problem in PROBLEMS:
        for method in METHODS:
            cfl = args.cfl if args.cfl is not None else METHOD_CFL[method]
            previous: dict[str, float | None] = {
                "l1": None,
                "l2": None,
                "linf": None,
            }
            for nx in args.resolutions:
                output_path = output_dir / f"{problem}_{method}_n{nx}.csv"
                command = [
                    str(executable),
                    f"--backend={args.backend}",
                    f"--problem={problem}",
                    f"--method={method}",
                    f"--nx={nx}",
                    "--boundary=periodic",
                    f"--final_time={args.final_time}",
                    f"--cfl={cfl}",
                    f"--output={output_path.relative_to(root).as_posix()}",
                ]
                print(
                    f"[report2-validation] problem={problem} "
                    f"method={method} backend={args.backend} "
                    f"nx={nx} cfl={cfl}"
                )
                subprocess.run(command, cwd=root, check=True)

                x, numerical = load_density(output_path)
                dx = (x_max - x_min) / nx
                exact = [
                    exact_density(
                        problem,
                        coordinate,
                        dx,
                        args.final_time,
                        velocity,
                        x_min,
                        x_max,
                    )
                    for coordinate in x
                ]
                l1, l2, linf = error_norms(numerical, exact)
                order_l1 = observed_order(previous["l1"], l1)
                order_l2 = observed_order(previous["l2"], l2)
                order_linf = observed_order(previous["linf"], linf)
                summary_rows.append(
                    {
                        "problem": problem,
                        "method": method,
                        "backend": args.backend,
                        "nx": nx,
                        "cfl": cfl,
                        "final_time": args.final_time,
                        "l1_rho": l1,
                        "order_l1": "" if order_l1 is None else order_l1,
                        "l2_rho": l2,
                        "order_l2": "" if order_l2 is None else order_l2,
                        "linf_rho": linf,
                        "order_linf": "" if order_linf is None else order_linf,
                    }
                )
                previous = {"l1": l1, "l2": l2, "linf": linf}

    summary_path = output_dir / "report2_high_order_convergence_1d.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    check_asymptotic_orders(summary_rows)
    print(f"saved: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
