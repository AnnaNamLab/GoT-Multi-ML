#!/usr/bin/env python3
"""Complete an expected-sequence CSV using IronThrone config files.

If the input CSV is missing targets required by the provided config files, this
helper derives those target sequences from the config files and writes a new
completed CSV. The original user-supplied CSV is never modified.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path

from build_expected_sequence_csv import expected_sequence_from_config, target_name_from_path


class ExpectedSequenceCsvError(RuntimeError):
    pass


def load_expected_sequence_csv(path: Path) -> OrderedDict[str, str]:
    entries: OrderedDict[str, str] = OrderedDict()
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for lineno, row in enumerate(reader, start=1):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) < 2:
                raise ExpectedSequenceCsvError(
                    f"{path}: line {lineno} must have at least two columns: target,expected_sequence"
                )
            target = row[0].strip()
            sequence = row[1].strip()
            if not target:
                raise ExpectedSequenceCsvError(f"{path}: line {lineno} has an empty target name")
            if target in entries:
                raise ExpectedSequenceCsvError(f"{path}: duplicate target in expected-sequence CSV: {target}")
            entries[target] = sequence
    return entries


def complete_expected_sequence_csv(
    config_files: list[Path],
    input_csv: Path,
    output_csv: Path,
    prefix: str,
    suffix: str,
) -> list[tuple[str, str]]:
    completed_entries = load_expected_sequence_csv(input_csv)
    missing_targets: list[tuple[str, str]] = []

    for config_path in config_files:
        target = target_name_from_path(config_path, prefix=prefix, suffix=suffix)
        had_target = target in completed_entries
        existing_sequence = completed_entries.get(target, "").strip()
        if existing_sequence:
            continue

        generated_sequence = expected_sequence_from_config(config_path)
        completed_entries[target] = generated_sequence
        reason = "blank-sequence" if had_target else "missing-target"
        missing_targets.append((target, reason))

    if not missing_targets:
        if output_csv.exists():
            output_csv.unlink()
        return []

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for target, sequence in completed_entries.items():
            writer.writerow([target, sequence])

    return missing_targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fill missing targets in an expected-sequence CSV using the provided config files. "
            "If no targets are missing, no output file is created."
        )
    )
    parser.add_argument("config_files", nargs="+", type=Path, help="Config files to check")
    parser.add_argument("--input-csv", required=True, type=Path, help="User-provided expected-sequence CSV")
    parser.add_argument("--output", required=True, type=Path, help="Path for the completed CSV if missing targets are found")
    parser.add_argument("--prefix", default="", help="Optional filename prefix to strip when deriving target names")
    parser.add_argument("--suffix", default=".config", help="Optional filename suffix to strip when deriving target names")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_csv = args.input_csv.resolve()
    output_csv = args.output.resolve()
    config_files = [path.resolve() for path in args.config_files]

    if not input_csv.is_file():
        raise ExpectedSequenceCsvError(f"Expected-sequence CSV not found: {input_csv}")

    missing_targets = complete_expected_sequence_csv(
        config_files=config_files,
        input_csv=input_csv,
        output_csv=output_csv,
        prefix=args.prefix,
        suffix=args.suffix,
    )

    if not missing_targets:
        print(f"Expected-sequence CSV already covers all requested targets. Using original file: {input_csv}")
        return 0

    for target, reason in missing_targets:
        if reason == "blank-sequence":
            print(
                f"Warning: target {target} has a blank expected sequence in {input_csv}. "
                f"Deriving it from the config file and writing a completed CSV."
            )
        else:
            print(
                f"Warning: target {target} was not found in {input_csv}. "
                f"Deriving it from the config file and writing a completed CSV."
            )

    print(f"Completed expected-sequence CSV written to: {output_csv}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ExpectedSequenceCsvError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
