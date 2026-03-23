#!/usr/bin/env python3
"""Build expected-sequence CSV files from linear IronThrone config files.

Rules implemented here:
- A linear config is expected to contain 4 rows:
  1. primer.sequence  start  end
  2. shared.sequence  start  end
  3. WT.sequence      start  end
  4. MUT.sequence     start  end
- The helper looks for the longest consecutive sequence block shared by WT and
  MUT based on positions, not fixed config-row order.
- If no shared consecutive block is longer than the shared sequence on line 2,
  it falls back to line 2.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True)
class ConfigRow:
    sequence: str
    start: int
    end: int


def parse_config_row(line: str, path: Path, lineno: int) -> ConfigRow:
    parts = line.rstrip("\n\r").split("\t")
    if len(parts) < 3:
        parts = line.rstrip("\n\r").split()
    if len(parts) < 3:
        raise ValueError(f"{path}: line {lineno} must have sequence, start, and end columns")
    sequence, start, end = parts[:3]
    return ConfigRow(sequence=sequence, start=int(start), end=int(end))


def load_linear_config(path: Path) -> list[ConfigRow]:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 4:
        raise ValueError(f"{path}: expected at least 4 non-empty lines for a linear config")
    return [parse_config_row(lines[i], path, i + 1) for i in range(4)]


def common_consecutive_rows(rows: Sequence[ConfigRow]) -> list[ConfigRow]:
    shared_rows: list[ConfigRow] = [rows[0], rows[1]]
    if (
        rows[2].sequence == rows[3].sequence
        and rows[2].start == rows[3].start
        and rows[2].end == rows[3].end
    ):
        shared_rows.append(rows[2])
    return sorted(shared_rows, key=lambda row: (row.start, row.end))


def joined_sequence(rows: Sequence[ConfigRow]) -> str:
    return "".join(row.sequence for row in rows)


def longest_common_consecutive_sequence(rows: Sequence[ConfigRow]) -> str:
    shared_rows = common_consecutive_rows(rows)
    best: list[ConfigRow] = []
    current: list[ConfigRow] = []

    for row in shared_rows:
        if not current or current[-1].end + 1 == row.start:
            current.append(row)
        else:
            if len(joined_sequence(current)) > len(joined_sequence(best)):
                best = current[:]
            current = [row]

    if len(joined_sequence(current)) > len(joined_sequence(best)):
        best = current[:]

    best_sequence = joined_sequence(best)
    shared_sequence = rows[1].sequence
    if len(best_sequence) > len(shared_sequence):
        return best_sequence
    return shared_sequence


def expected_sequence_from_rows(rows: Sequence[ConfigRow]) -> str:
    return longest_common_consecutive_sequence(rows)


def expected_sequence_from_config(path: Path) -> str:
    return expected_sequence_from_rows(load_linear_config(path))


def target_name_from_path(path: Path, prefix: str = "", suffix: str = ".config") -> str:
    name = path.name
    if suffix and name.endswith(suffix):
        name = name[: -len(suffix)]
    if prefix and name.startswith(prefix):
        name = name[len(prefix) :]
    return name


def write_expected_sequence_csv(
    config_files: Iterable[Path],
    output: Path,
    prefix: str = "",
    suffix: str = ".config",
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for config_path in config_files:
            writer.writerow(
                [
                    target_name_from_path(config_path, prefix=prefix, suffix=suffix),
                    expected_sequence_from_config(config_path),
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build expected_sequence.csv from linear IronThrone config files."
    )
    parser.add_argument("config_files", nargs="+", type=Path, help="Config files to process")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output CSV path. If omitted, CSV is written to stdout.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        default="",
        help="Optional filename prefix to strip when deriving target names.",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        default=".config",
        help="Optional filename suffix to strip when deriving target names.",
    )
    parser.add_argument(
        "--sequence-only",
        action="store_true",
        help="Print only the expected sequence for a single config file.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config_files = [path.resolve() for path in args.config_files]

    if args.sequence_only:
        if len(config_files) != 1:
            raise SystemExit("--sequence-only requires exactly one config file")
        print(expected_sequence_from_config(config_files[0]))
        return 0

    if args.output is not None:
        write_expected_sequence_csv(
            config_files,
            output=args.output.resolve(),
            prefix=args.prefix,
            suffix=args.suffix,
        )
        return 0

    writer = csv.writer(sys.stdout, lineterminator="\n")
    for config_path in config_files:
        writer.writerow(
            [
                target_name_from_path(config_path, prefix=args.prefix, suffix=args.suffix),
                expected_sequence_from_config(config_path),
            ]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
