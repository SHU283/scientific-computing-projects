#!/usr/bin/env python3
"""Audit CFL metadata in canonical Report 2 HLLC result tables."""

from __future__ import annotations

import csv
import ast
from pathlib import Path


SMOOTH_CFL = {"weno2": 0.4, "weno3": 0.4, "ader3": 0.25}
DISCONTINUOUS_CFL = {"weno2": 0.2, "weno3": 0.2, "ader3": 0.14}
SKIP_PARTS = {
    "archive",
    "archive_old_cuda_checks_before_weno2_fix",
    "analysis_candidates",
    "ader3_cfl_probe",
    "ader3_cfl012_confirm",
    "ader3_diagnostics",
    "ader3_discontinuous_cfl_probe",
    "ader3_quadrant_cfl_probe",
    "backend_consistency_2d_lovelace_repeat",
}

SCRIPT_MAPS = {
    "audit_backend_alignment_2d.py": {
        "SMOOTH_CFL": SMOOTH_CFL,
        "DISCONTINUOUS_CFL": DISCONTINUOUS_CFL,
    },
    "run_report2_backend_consistency_2d.py": {
        "SMOOTH_CFL": SMOOTH_CFL,
        "DISCONTINUOUS_CFL": DISCONTINUOUS_CFL,
    },
    "run_report2_discontinuous_runtime_2d.py": {"METHOD_CFL": DISCONTINUOUS_CFL},
    "run_report2_gpu_limit_extension_2d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_report2_high_resolution_fields_1d.py": {
        "SMOOTH_CFL": SMOOTH_CFL,
        "DISCONTINUOUS_CFL": DISCONTINUOUS_CFL,
    },
    "run_report2_local_smoke.py": {
        "SMOOTH_CFL": SMOOTH_CFL,
        "DISCONTINUOUS_CFL": DISCONTINUOUS_CFL,
    },
    "run_report2_openmp_thread_scaling_2d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_report2_phase_profile_2d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_report2_runtime_scaling_2d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_report2_validation_1d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_report2_validation_2d.py": {"METHOD_CFL": SMOOTH_CFL},
    "run_stageB_high_order_1d_checks.py": {"SMOOTH_CFL": SMOOTH_CFL},
    "run_high_order_runtime_1d.py": {"REPORT2_CFL": SMOOTH_CFL},
    "run_high_order_single_cases_1d.py": {"REPORT2_CFL": SMOOTH_CFL},
}


def expected_cfl(problem: str, method: str) -> float | None:
    if method not in SMOOTH_CFL:
        return None
    if problem in {"smooth_advection", "gaussian_advection", "gaussian_advection_2d"}:
        return SMOOTH_CFL[method]
    if problem in {"sod", "quadrant_riemann_2d", "shock_bubble_2d"}:
        return DISCONTINUOUS_CFL[method]
    return None


def assigned_literal(path: Path, variable: str) -> object | None:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == variable for target in node.targets):
            return ast.literal_eval(node.value)
    return None


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    package_root = Path(__file__).resolve().parents[3]
    result_root = package_root / "results" / "hllc"
    failures: list[str] = []
    checked = 0
    checked_maps = 0

    script_root = root / "scripts"
    for filename, expected_maps in SCRIPT_MAPS.items():
        path = script_root / filename
        for variable, expected in expected_maps.items():
            checked_maps += 1
            actual = assigned_literal(path, variable)
            if actual != expected:
                failures.append(
                    f"scripts/{filename}: {variable}={actual!r}, expected={expected!r}"
                )

    for path in sorted(result_root.rglob("*.csv")):
        if any(part in SKIP_PARTS for part in path.parts):
            continue
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or not {"problem", "method", "cfl"}.issubset(
                reader.fieldnames
            ):
                continue
            for line_number, row in enumerate(reader, start=2):
                expected = expected_cfl(row["problem"], row["method"])
                if expected is None:
                    continue
                checked += 1
                actual = float(row["cfl"])
                if abs(actual - expected) > 1.0e-12:
                    failures.append(
                        f"{path.relative_to(package_root)}:{line_number}: "
                        f"{row['problem']}/{row['method']} CFL={actual}, expected={expected}"
                    )

    print(f"checked active CFL maps: {checked_maps}")
    print(f"checked canonical CFL rows: {checked}")
    if failures:
        print("CFL POLICY FAIL")
        print("\n".join(failures))
        return 1
    print("CFL POLICY PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
