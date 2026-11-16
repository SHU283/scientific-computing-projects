#!/usr/bin/env python3
"""Merge the corrected WENO-2 OpenMP refresh into canonical Report 2 tables.

Only rows whose experiment keys occur in the refresh files are replaced.  The
original canonical CSV files are copied to a dated archive before any write.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime
from pathlib import Path


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise RuntimeError(f"CSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def key_for(row: dict[str, str], fields: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(row[field] for field in fields)


def merge_one(
    canonical: Path,
    refresh: Path,
    key_fields: tuple[str, ...],
    archive_dir: Path,
    allowed_resolutions: set[int] | None = None,
) -> int:
    fields, rows = read_csv(canonical)
    refresh_fields, refresh_rows = read_csv(refresh)
    if fields != refresh_fields:
        raise RuntimeError(
            f"Header mismatch between {canonical.name} and {refresh.name}"
        )

    if allowed_resolutions is not None:
        refresh_rows = [
            row for row in refresh_rows if int(row["nx"]) in allowed_resolutions
        ]
    replacements = {key_for(row, key_fields): row for row in refresh_rows}
    if len(replacements) != len(refresh_rows):
        raise RuntimeError(f"Duplicate refresh keys in {refresh}")

    canonical_keys = {key_for(row, key_fields) for row in rows}
    missing = set(replacements) - canonical_keys
    if missing:
        raise RuntimeError(f"Refresh keys are absent from {canonical}: {sorted(missing)}")

    merged = [replacements.get(key_for(row, key_fields), row) for row in rows]
    changed = sum(
        old != new for old, new in zip(rows, merged)
    )
    if changed != len(replacements):
        raise RuntimeError(
            f"Expected {len(replacements)} changed rows in {canonical}, got {changed}"
        )

    archive_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(canonical, archive_dir / canonical.name)
    write_csv(canonical, fields, merged)
    print(f"updated {canonical}: {changed} row(s)")
    return changed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="report2-hllc repository root",
    )
    parser.add_argument(
        "--archive-label",
        default=datetime.now().strftime("pre_weno2_openmp_refresh_%Y%m%d_%H%M%S"),
    )
    args = parser.parse_args()

    performance = args.root / "results" / "report2_hllc" / "performance"
    refresh_dir = performance / "weno2_openmp_refresh"
    archive_dir = performance / "archive" / args.archive_label

    smooth = performance / "runtime_accuracy_2d_report2_hllc_lovelace_a30.csv"
    discontinuous = performance / "discontinuous"
    thread = performance / "openmp_thread_scaling"
    phase = performance / "phase_profile"

    total = 0
    total += merge_one(
        smooth,
        refresh_dir / "smooth_weno2_openmp_refresh.csv",
        ("problem", "method", "backend", "nx", "ny"),
        archive_dir,
    )
    for target, resolutions in (
        (
            discontinuous
            / "runtime_discontinuous_2d_report2_hllc_lovelace_a30_50_100_200.csv",
            {50, 100, 200},
        ),
        (
            discontinuous
            / "runtime_discontinuous_2d_report2_hllc_lovelace_a30_n600.csv",
            {600},
        ),
        (
            discontinuous
            / "runtime_discontinuous_2d_report2_hllc_lovelace_a30_50_100_200_600.csv",
            {50, 100, 200, 600},
        ),
    ):
        total += merge_one(
            target,
            refresh_dir / "discontinuous_weno2_openmp_refresh.csv",
            ("problem", "method", "backend", "nx", "ny"),
            archive_dir,
            resolutions,
        )
    total += merge_one(
        thread / "openmp_thread_scaling_report2_hllc_lovelace_a30.csv",
        refresh_dir / "openmp_thread_scaling_weno2_refresh.csv",
        ("problem", "method", "backend", "nx", "ny", "threads"),
        archive_dir,
    )
    total += merge_one(
        phase / "phase_profile_2d_hllc_lovelace_a30.csv",
        refresh_dir / "phase_profile_weno2_openmp_refresh.csv",
        ("problem", "method", "backend", "nx", "ny", "repeat"),
        archive_dir,
    )
    print(f"archive: {archive_dir}")
    print(f"total replaced rows: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
