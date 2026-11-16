#!/usr/bin/env python3
"""Verify the standalone report's figure assets and formal CFL statements."""

from __future__ import annotations

import filecmp
import re
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[3]
REPORT = PACKAGE / "report" / "standalone_report2"
REPORT_FIGURES = REPORT / "figures"
CURATED_FIGURES = PACKAGE / "figures" / "formal_report"


def main() -> int:
    failures: list[str] = []
    included: set[str] = set()
    pattern = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")

    if not REPORT.is_dir():
        print(f"FORMAL REPORT ASSET AUDIT FAIL\nmissing report directory: {REPORT}")
        return 1

    for tex in REPORT.rglob("*.tex"):
        included.update(pattern.findall(tex.read_text(encoding="utf-8-sig")))

    checked = 0
    for name in sorted(included):
        report_asset = REPORT_FIGURES / name
        curated_asset = CURATED_FIGURES / name
        if not report_asset.is_file():
            failures.append(f"missing standalone report asset: {name}")
            continue
        if not curated_asset.is_file():
            failures.append(f"missing curated formal asset: {name}")
            continue
        checked += 1
        if not filecmp.cmp(report_asset, curated_asset, shallow=False):
            failures.append(f"content mismatch between report and curated asset: {name}")

    unreferenced = sorted(
        path.name
        for path in REPORT_FIGURES.iterdir()
        if path.is_file() and path.name != "University_Crest.pdf" and path.name not in included
    )
    if unreferenced:
        failures.append(f"unreferenced standalone report figures: {unreferenced}")

    methodology = (REPORT / "sections" / "03_numerical_methods.tex").read_text(
        encoding="utf-8-sig"
    )
    required_policy = (
        "Smooth runs use CFL 0.4 for WENO-2/WENO-3 and 0.25 for ADER-3.",
        "Discontinuous WENO runs use 0.20.",
        "common ADER-3 value 0.14",
    )
    for statement in required_policy:
        if statement not in methodology:
            failures.append(f"formal CFL statement missing: {statement}")

    print(f"formal assets referenced: {len(included)}")
    print(f"report/curated file pairs checked: {checked}")
    if failures:
        print("FORMAL REPORT ASSET AUDIT FAIL")
        print("\n".join(failures))
        return 1
    print("FORMAL REPORT ASSET AUDIT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
