#!/usr/bin/env python3
"""Recheck WENO-2/WENO-3 smooth 2D convergence on one backend."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
from pathlib import Path


METHODS = ("weno2", "weno3")
DEFAULT_RESOLUTIONS = (25, 50, 100, 200, 400)
ERROR_PATTERN = re.compile(
    r"2D density error L1\s+:\s+(?P<l1>\S+).*?"
    r"2D density error L2\s+:\s+(?P<l2>\S+).*?"
    r"2D density error Linf\s+:\s+(?P<linf>\S+)",
    re.DOTALL,
)
STEP_PATTERN = re.compile(r"time steps\s+:\s+(?P<steps>\d+)")


def observed_order(coarse: float | None, fine: float) -> float | None:
    if coarse is None or coarse <= 0.0 or fine <= 0.0:
        return None
    return math.log(coarse / fine, 2.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", default="build-report2-hllc-cuda")
    parser.add_argument(
        "--output-dir",
        default="results/report2_hllc/validation_2d_weno_cuda_recheck",
    )
    parser.add_argument(
        "--resolutions", nargs="+", type=int, default=DEFAULT_RESOLUTIONS
    )
    parser.add_argument(
        "--methods", nargs="+", choices=METHODS, default=METHODS
    )
    parser.add_argument(
        "--backend", choices=("serial", "openmp", "cuda"), default="cuda"
    )
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.4)
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    executable = root / args.build_dir / "final_project"
    output_dir = root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    if not executable.exists():
        raise SystemExit(f"solver executable not found: {executable}")

    rows: list[dict[str, str | int | float]] = []
    for method in args.methods:
        previous: dict[str, float | None] = {"l1": None, "l2": None, "linf": None}
        for n in args.resolutions:
            command = [
                str(executable),
                f"--backend={args.backend}",
                "--problem=gaussian_advection_2d",
                f"--method={method}",
                f"--nx={n}",
                f"--ny={n}",
                "--boundary=periodic",
                f"--final_time={args.final_time}",
                f"--cfl={args.cfl}",
                "--output=none",
            ]
            print(
                f"[weno-convergence-recheck] method={method} "
                f"backend={args.backend} n={n} cfl={args.cfl}",
                flush=True,
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
            step_match = STEP_PATTERN.search(completed.stdout)
            values = {
                "l1": float(match.group("l1")),
                "l2": float(match.group("l2")),
                "linf": float(match.group("linf")),
            }
            row: dict[str, str | int | float] = {
                "method": method,
                "backend": args.backend,
                "n": n,
                "cfl": args.cfl,
                "final_time": args.final_time,
                "steps": int(step_match.group("steps")) if step_match else "",
                "l1_rho": values["l1"],
                "order_l1": "",
                "l2_rho": values["l2"],
                "order_l2": "",
                "linf_rho": values["linf"],
                "order_linf": "",
            }
            for name, value in values.items():
                order = observed_order(previous[name], value)
                row[f"order_{name}"] = "" if order is None else order
                previous[name] = value
            rows.append(row)
            print(
                f"  L1={values['l1']:.8e} "
                f"order={row['order_l1'] if row['order_l1'] != '' else '-'}",
                flush=True,
            )

    summary_path = output_dir / "weno2_weno3_convergence_cuda_recheck.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    expected = {"weno2": (1.8, 2.2), "weno3": (2.7, 3.3)}
    for method, bounds in expected.items():
        method_rows = [row for row in rows if row["method"] == method]
        if not method_rows:
            continue
        final_order = float(method_rows[-1]["order_l1"])
        if not bounds[0] <= final_order <= bounds[1]:
            raise RuntimeError(
                f"{method} finest-grid L1 order outside expected range: "
                f"{final_order:.6f}"
            )
        print(
            f"{method} finest-grid L1 order: PASS ({final_order:.6f})",
            flush=True,
        )

    print(f"saved: {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
