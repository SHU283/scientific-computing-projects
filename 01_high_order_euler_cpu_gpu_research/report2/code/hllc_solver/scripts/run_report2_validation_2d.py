#!/usr/bin/env python3
"""Run isolated 2D convergence checks for the Report 2 methods."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
DEFAULT_RESOLUTIONS = (25, 50, 100, 200, 400)
METHOD_CFL = {
    "weno2": 0.4,
    "weno3": 0.4,
    "ader3": 0.25,
}
ERROR_PATTERN = re.compile(
    r"2D density error L1\s+:\s+(?P<l1>\S+).*?"
    r"2D density error L2\s+:\s+(?P<l2>\S+).*?"
    r"2D density error Linf\s+:\s+(?P<linf>\S+)",
    re.DOTALL,
)


def observed_order(coarse: float | None, fine: float) -> float | None:
    if coarse is None or coarse <= 0.0 or fine <= 0.0:
        return None
    return math.log(coarse / fine, 2.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2")
    parser.add_argument("--output-dir", default="results/report2/validation_2d")
    parser.add_argument(
        "--resolutions", nargs="+", type=int, default=DEFAULT_RESOLUTIONS
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=METHODS,
        default=METHODS,
        help="Methods to validate (default: all Report 2 methods).",
    )
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument(
        "--backend", choices=("serial", "openmp", "cuda"), default="serial"
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

    rows: list[dict[str, str | int | float]] = []
    for method in args.methods:
        cfl = args.cfl if args.cfl is not None else METHOD_CFL[method]
        previous: dict[str, float | None] = {"l1": None, "l2": None, "linf": None}
        for n in args.resolutions:
            output_path = output_dir / f"gaussian_advection_2d_{method}_n{n}.csv"
            command = [
                str(executable),
                f"--backend={args.backend}",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={n}",
                f"--ny={n}",
                "--boundary=periodic",
                f"--final_time={args.final_time}",
                f"--cfl={cfl}",
                f"--output={output_path.relative_to(root).as_posix()}",
            ]
            print(
                f"[report2-validation-2d] method={method} "
                f"backend={args.backend} n={n} cfl={cfl}"
            )
            completed = subprocess.run(
                command, cwd=root, check=True, capture_output=True, text=True
            )
            match = ERROR_PATTERN.search(completed.stdout)
            if not match:
                raise RuntimeError(
                    f"could not read error norms for {method} n={n}\n"
                    f"{completed.stdout}\n{completed.stderr}"
                )

            l1 = float(match.group("l1"))
            l2 = float(match.group("l2"))
            linf = float(match.group("linf"))
            row: dict[str, str | int | float] = {
                "method": method,
                "backend": args.backend,
                "n": n,
                "cfl": cfl,
                "final_time": args.final_time,
                "l1_rho": l1,
                "order_l1": "",
                "l2_rho": l2,
                "order_l2": "",
                "linf_rho": linf,
                "order_linf": "",
            }
            for name, value in (("l1", l1), ("l2", l2), ("linf", linf)):
                order = observed_order(previous[name], value)
                row[f"order_{name}"] = "" if order is None else order
                previous[name] = value
            rows.append(row)

    summary_path = output_dir / "report2_high_order_convergence_2d.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    expected_orders = {
        "weno2": (1.8, 2.2),
        "weno3": (2.7, 3.3),
        "ader3": (2.7, 3.3),
    }
    for method, bounds in expected_orders.items():
        if method not in args.methods:
            continue
        method_rows = [row for row in rows if row["method"] == method]
        final_order = float(method_rows[-1]["order_l1"])
        if not bounds[0] <= final_order <= bounds[1]:
            raise RuntimeError(
                f"2D {method} L1 convergence check failed: order={final_order:.3f}"
            )
        print(
            f"2D {method} convergence check: PASS "
            f"(L1 order={final_order:.3f})"
        )
    print(f"saved: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
