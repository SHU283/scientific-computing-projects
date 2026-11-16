#!/usr/bin/env python3
"""Run a compact local regression suite for the Report 2 solvers."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path


FIELDS_1D = ("rho", "u", "p", "energy")
FIELDS_2D = ("rho", "u", "v", "p", "energy")
SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
DISCONTINUOUS_CFL = {"weno2": 0.2, "weno3": 0.2, "ader3": 0.14}


def run_case(executable: Path, root: Path, arguments: list[str]) -> None:
    command = [str(executable), *arguments]
    print(f"\n$ {' '.join(command)}")
    subprocess.run(command, cwd=root, check=True)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def check_admissible(path: Path, fields: tuple[str, ...]) -> None:
    rows = load_rows(path)
    if not rows:
        raise RuntimeError(f"empty result file: {path}")
    for row_index, row in enumerate(rows):
        values = [float(row[field]) for field in fields]
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"non-finite state in {path} at row {row_index}")
        if float(row["rho"]) <= 0.0 or float(row["p"]) <= 0.0:
            raise RuntimeError(f"non-positive state in {path} at row {row_index}")


def compare_results(
    first: Path,
    second: Path,
    fields: tuple[str, ...],
    tolerance: float,
) -> dict[str, float]:
    first_rows = load_rows(first)
    second_rows = load_rows(second)
    if len(first_rows) != len(second_rows):
        raise RuntimeError(
            f"row-count mismatch: {first.name}={len(first_rows)}, "
            f"{second.name}={len(second_rows)}"
        )

    differences = {
        field: max(
            abs(float(left[field]) - float(right[field]))
            for left, right in zip(first_rows, second_rows)
        )
        for field in fields
    }
    failed = {field: value for field, value in differences.items() if value > tolerance}
    if failed:
        details = ", ".join(f"{field}={value:.3e}" for field, value in failed.items())
        raise RuntimeError(
            f"{first.name} and {second.name} differ above {tolerance:.1e}: {details}"
        )
    return differences


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2")
    parser.add_argument("--output-dir", default="results/report2/local_smoke")
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    executable = root / args.build_dir / "final_project"
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    comparisons: list[tuple[str, dict[str, float]]] = []

    for method in ("weno3", "ader3"):
        serial = output_dir / f"smooth_{method}_serial_n120.csv"
        openmp = output_dir / f"smooth_{method}_openmp_n120.csv"
        common = [
            "--problem=smooth_advection",
            f"--method={method}",
            "--nx=120",
            "--boundary=periodic",
            "--final_time=0.05",
            f"--cfl={SMOOTH_CFL[method]}",
        ]
        run_case(
            executable,
            root,
            ["--backend=serial", *common, f"--output={serial.relative_to(root).as_posix()}"],
        )
        run_case(
            executable,
            root,
            ["--backend=openmp", *common, f"--output={openmp.relative_to(root).as_posix()}"],
        )
        check_admissible(serial, FIELDS_1D)
        check_admissible(openmp, FIELDS_1D)
        comparisons.append(
            (
                f"1D {method} serial/OpenMP",
                compare_results(serial, openmp, FIELDS_1D, args.tolerance),
            )
        )

    for method in ("weno3", "ader3"):
        serial_2d = output_dir / f"gaussian_{method}_serial_n40.csv"
        openmp_2d = output_dir / f"gaussian_{method}_openmp_n40.csv"
        common_2d = [
            "--problem=gaussian_advection_2d",
            f"--method={method}",
            "--nx=40",
            "--ny=40",
            "--boundary=periodic",
            "--final_time=0.02",
            f"--cfl={SMOOTH_CFL[method]}",
        ]
        run_case(
            executable,
            root,
            [
                "--backend=serial",
                *common_2d,
                f"--output={serial_2d.relative_to(root).as_posix()}",
            ],
        )
        run_case(
            executable,
            root,
            [
                "--backend=openmp",
                *common_2d,
                f"--output={openmp_2d.relative_to(root).as_posix()}",
            ],
        )
        check_admissible(serial_2d, FIELDS_2D)
        check_admissible(openmp_2d, FIELDS_2D)
        comparisons.append(
            (
                f"2D {method} serial/OpenMP",
                compare_results(serial_2d, openmp_2d, FIELDS_2D, args.tolerance),
            )
        )

    discontinuous_cases = (
        (
            "sod",
            "ader3",
            200,
            None,
            0.08,
            DISCONTINUOUS_CFL["ader3"],
            output_dir / "sod_ader3_n200.csv",
        ),
        (
            "quadrant_riemann_2d",
            "weno3",
            64,
            64,
            0.01,
            0.2,
            output_dir / "quadrant_weno3_n64.csv",
        ),
        (
            "shock_bubble_2d",
            "weno3",
            64,
            64,
            0.01,
            0.2,
            output_dir / "shock_bubble_weno3_n64.csv",
        ),
        (
            "quadrant_riemann_2d",
            "ader3",
            40,
            40,
            0.01,
            0.14,
            output_dir / "quadrant_ader3_n40.csv",
        ),
        (
            "shock_bubble_2d",
            "ader3",
            40,
            40,
            0.01,
            0.14,
            output_dir / "shock_bubble_ader3_n40.csv",
        ),
    )
    for problem, method, nx, ny, final_time, cfl, output in discontinuous_cases:
        arguments = [
            "--backend=serial",
            f"--problem={problem}",
            f"--method={method}",
            f"--nx={nx}",
            "--boundary=outflow",
            f"--final_time={final_time}",
            f"--cfl={cfl}",
            f"--output={output.relative_to(root).as_posix()}",
        ]
        if ny is not None:
            arguments.append(f"--ny={ny}")
        run_case(executable, root, arguments)
        check_admissible(output, FIELDS_2D if ny is not None else FIELDS_1D)

    summary = output_dir / "local_smoke_summary.txt"
    lines = ["Report 2 local smoke tests: PASS"]
    for label, differences in comparisons:
        lines.append(label)
        lines.extend(
            f"  max_abs_{field}={value:.17e}"
            for field, value in differences.items()
        )
    summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n{summary.read_text(encoding='utf-8')}", end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
