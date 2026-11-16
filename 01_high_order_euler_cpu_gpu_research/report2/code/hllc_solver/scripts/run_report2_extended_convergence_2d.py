#!/usr/bin/env python3
"""Run the fine-grid 2D HLLC convergence extension without field output."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path


METHODS = ("weno2", "weno3", "ader3")
METHOD_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
EXPECTED_FINAL_ORDER = {
    "weno2": (1.8, 2.2),
    "weno3": (2.7, 3.3),
    "ader3": (2.7, 3.3),
}
ERROR_PATTERN = re.compile(
    r"2D density error L1\s+:\s+(?P<l1>\S+).*?"
    r"2D density error L2\s+:\s+(?P<l2>\S+).*?"
    r"2D density error Linf\s+:\s+(?P<linf>\S+)",
    re.DOTALL,
)


def observed_order(
    coarse_n: int | None,
    coarse_error: float | None,
    fine_n: int,
    fine_error: float,
) -> float | None:
    if (
        coarse_n is None
        or coarse_error is None
        or coarse_error <= 0.0
        or fine_error <= 0.0
    ):
        return None
    return math.log(coarse_error / fine_error) / math.log(fine_n / coarse_n)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2-hllc-extended")
    parser.add_argument(
        "--output-dir",
        default=(
            "results/report2_hllc/extended_convergence_1600/"
            "validation_2d_cuda"
        ),
    )
    parser.add_argument(
        "--resolutions", nargs="+", type=int, default=(400, 800, 1600)
    )
    parser.add_argument("--final-time", type=float, default=0.1)
    args = parser.parse_args()

    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2:
        raise SystemExit("at least two resolutions are required")

    root = Path(__file__).resolve().parent.parent
    executable = root / args.build_dir / "final_project"
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    rows: list[dict[str, str | int | float]] = []
    for method in METHODS:
        cfl = METHOD_CFL[method]
        previous_n: int | None = None
        previous = {"l1": None, "l2": None, "linf": None}
        for n in resolutions:
            command = [
                str(executable),
                "--backend=cuda",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={n}",
                f"--ny={n}",
                "--boundary=periodic",
                f"--final_time={args.final_time}",
                f"--cfl={cfl}",
                "--output=none",
            ]
            print(f"[extended-2d] method={method} n={n} cfl={cfl}", flush=True)
            completed = subprocess.run(
                command,
                cwd=root,
                check=True,
                capture_output=True,
                text=True,
            )
            match = ERROR_PATTERN.search(completed.stdout)
            if not match:
                raise RuntimeError(
                    f"could not parse error norms for {method} n={n}\n"
                    f"{completed.stdout}\n{completed.stderr}"
                )

            values = {
                "l1": float(match.group("l1")),
                "l2": float(match.group("l2")),
                "linf": float(match.group("linf")),
            }
            orders = {
                name: observed_order(previous_n, previous[name], n, value)
                for name, value in values.items()
            }
            rows.append(
                {
                    "method": method,
                    "backend": "cuda",
                    "n": n,
                    "cfl": cfl,
                    "final_time": args.final_time,
                    "l1_rho": values["l1"],
                    "order_l1": "" if orders["l1"] is None else orders["l1"],
                    "l2_rho": values["l2"],
                    "order_l2": "" if orders["l2"] is None else orders["l2"],
                    "linf_rho": values["linf"],
                    "order_linf": (
                        "" if orders["linf"] is None else orders["linf"]
                    ),
                }
            )
            previous_n = n
            previous = values
            order_text = "-" if orders["l1"] is None else f"{orders['l1']:.6f}"
            print(
                f"  L1={values['l1']:.8e} order={order_text}",
                flush=True,
            )

    output_path = output_dir / "extended_convergence_2d_cuda_400_800_1600.csv"
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    failures: list[str] = []
    for method in METHODS:
        final_row = [row for row in rows if row["method"] == method][-1]
        final_order = float(final_row["order_l1"])
        lower, upper = EXPECTED_FINAL_ORDER[method]
        if not lower <= final_order <= upper:
            failures.append(
                f"{method}: final L1 order {final_order:.4f} outside "
                f"[{lower:.1f}, {upper:.1f}]"
            )

    print(f"saved: {output_path}")
    if failures:
        raise RuntimeError("fine-grid checks failed:\n  " + "\n  ".join(failures))
    print("2D 800-to-1600 fine-grid checks: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
