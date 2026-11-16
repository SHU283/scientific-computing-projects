#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMPARE = ROOT / "scripts" / "compare_csv_results.py"


CASES = [
    {
        "problem": "smooth_advection",
        "nx": 200,
        "boundary": "periodic",
        "final_time": 0.3,
    },
    {
        "problem": "gaussian_advection",
        "nx": 200,
        "boundary": "periodic",
        "final_time": 0.3,
    },
    {
        "problem": "sod",
        "nx": 200,
        "boundary": "outflow",
        "final_time": 0.2,
    },
    {
        "problem": "gaussian_advection_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "periodic",
        "final_time": 0.1,
    },
    {
        "problem": "quadrant_riemann_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "outflow",
        "final_time": 0.2,
    },
    {
        "problem": "shock_bubble_2d",
        "nx": 100,
        "ny": 100,
        "boundary": "outflow",
        "final_time": 0.05,
    },
]


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    print("$", " ".join(cmd))
    completed = subprocess.run(cmd, cwd=ROOT, check=True, text=True, capture_output=True)
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")
    print()
    return completed


def find_executable() -> Path:
    candidates = [
        ROOT / "build-cuda" / "final_project",
        ROOT / "build" / "final_project",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit("No final_project executable found in build/ or build-cuda/.")


def solver_cmd(executable: Path, backend: str, case: dict[str, object], output: Path) -> list[str]:
    cmd = [
        str(executable),
        f"--backend={backend}",
        f"--problem={case['problem']}",
        "--method=weno2",
        f"--nx={case['nx']}",
        f"--boundary={case['boundary']}",
        f"--final_time={case['final_time']}",
        f"--output={output.as_posix()}",
    ]
    if "ny" in case:
        cmd.insert(5, f"--ny={case['ny']}")
    return cmd


def main() -> int:
    executable = find_executable()
    if "build-cuda" not in executable.as_posix():
        raise SystemExit(
            "Found only a non-CUDA executable. Please build with `cmake -S . -B build-cuda -DFP_ENABLE_CUDA=ON` first."
        )

    stage_b = ROOT / "results" / "stageB"
    stage_b.mkdir(parents=True, exist_ok=True)

    for case in CASES:
        problem = str(case["problem"])
        serial_csv = stage_b / f"{problem}_serial_weno2.csv"
        cuda_csv = stage_b / f"{problem}_cuda_weno2.csv"

        run(solver_cmd(executable, "serial", case, serial_csv))
        completed = run(solver_cmd(executable, "cuda", case, cuda_csv))
        if "CUDA backend not compiled into this build." in completed.stdout or not cuda_csv.exists():
            raise SystemExit(
                "The selected executable does not contain CUDA kernels. Rebuild with FP_ENABLE_CUDA=ON and rerun."
            )
        run(["python3", str(COMPARE), str(serial_csv), str(cuda_csv)])

    print("Report 1 WENO-2 CUDA backend checks finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

