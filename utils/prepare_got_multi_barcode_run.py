#!/usr/bin/env python3
"""Prepare the BARCODE second pass for GoT-Multi.

This helper keeps the mutation-side input preparation from the first IronThrone-Multi
run, then rewrites the cached processed_input params CSV so the second pass uses the
barcode configs target-by-target.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path


class ParamsError(RuntimeError):
    pass


def parse_params_file(path: Path) -> tuple[list[str], OrderedDict[str, list[str]]]:
    order: list[str] = []
    sections: OrderedDict[str, list[str]] = OrderedDict()
    current: str | None = None

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n").rstrip("\r")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                current = stripped[1:-1].strip()
                if current not in sections:
                    sections[current] = []
                    order.append(current)
                continue
            if current is None:
                raise ParamsError(f"Found content before any section header in {path}: {line}")
            sections[current].append(stripped)

    return order, sections


def parse_ironthrone_args(lines: list[str], source: Path) -> OrderedDict[str, str]:
    args: OrderedDict[str, str] = OrderedDict()
    for line in lines:
        parts = [part.strip() for part in line.split(",", 1)]
        if len(parts) != 2 or not parts[0]:
            raise ParamsError(f"Malformed [ironthrone-args] line in {source}: {line}")
        args[parts[0]] = parts[1]
    return args


def expand_config_entries(entries: list[str], source: Path) -> list[Path]:
    if not entries:
        raise ParamsError(f"[config-files] is empty in {source}")

    dirs: list[Path] = []
    files: list[Path] = []
    for entry in entries:
        path = Path(entry).expanduser()
        if path.is_dir():
            dirs.append(path)
        else:
            files.append(path)

    if len(dirs) > 1:
        raise ParamsError(f"[config-files] accepts at most one directory in {source}")
    if dirs and files:
        raise ParamsError(f"[config-files] cannot mix a directory with individual files in {source}")

    if dirs:
        config_dir = dirs[0]
        configs = sorted(path.resolve() for path in config_dir.iterdir() if path.is_file() and path.suffix == ".config")
        if not configs:
            raise ParamsError(f"No .config files found in directory: {config_dir}")
        return configs

    resolved_files = [path.resolve() for path in files]
    missing = [str(path) for path in resolved_files if not path.is_file()]
    if missing:
        raise ParamsError(f"Config file(s) not found in {source}: {', '.join(missing)}")
    return resolved_files


def extract_target_name(config_path: Path, prefixes: list[str], suffix: str) -> str:
    name = config_path.name
    if suffix and not name.endswith(suffix):
        raise ParamsError(
            f"Config file does not match the expected suffix {suffix!r}: {config_path}"
        )

    stem = name[:-len(suffix)] if suffix else config_path.stem
    for prefix in prefixes:
        if prefix and stem.startswith(prefix):
            return stem[len(prefix) :]
    return stem


def find_processed_params_csv(processed_input_dir: Path) -> Path:
    matches = sorted(processed_input_dir.glob("*ironethrone_params.csv"))
    if not matches:
        raise ParamsError(f"No *ironethrone_params.csv found in {processed_input_dir}")
    if len(matches) > 1:
        joined = ", ".join(str(path) for path in matches)
        raise ParamsError(f"Multiple *ironethrone_params.csv files found in {processed_input_dir}: {joined}")
    return matches[0]


def write_params_file(order: list[str], sections: OrderedDict[str, list[str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="\n") as handle:
        for idx, section in enumerate(order):
            handle.write(f"[{section}]\n")
            for line in sections[section]:
                handle.write(f"{line}\n")
            if idx != len(order) - 1:
                handle.write("\n")


def rewrite_processed_params_csv(csv_path: Path, barcode_config_by_target: dict[str, Path]) -> tuple[list[str], list[str]]:
    with csv_path.open(newline="") as handle:
        rows = list(csv.reader(handle))

    updated_rows: list[list[str]] = []
    seen_targets: list[str] = []
    missing_targets: list[str] = []

    for row in rows:
        if len(row) != 5:
            raise ParamsError(f"Expected 5 columns in {csv_path}, got {len(row)} columns: {row}")
        target = row[0].strip()
        seen_targets.append(target)
        if target not in barcode_config_by_target:
            missing_targets.append(target)
            continue
        updated_row = row[:]
        updated_row[4] = str(barcode_config_by_target[target])
        updated_rows.append(updated_row)

    if missing_targets:
        missing_str = ", ".join(sorted(set(missing_targets)))
        raise ParamsError(
            f"BARCODE configs are missing the following targets required by processed_input: {missing_str}"
        )

    extra_targets = sorted(set(barcode_config_by_target) - set(seen_targets))
    if extra_targets:
        extra_str = ", ".join(extra_targets)
        raise ParamsError(
            f"BARCODE configs contain targets not present in processed_input: {extra_str}"
        )

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(updated_rows)

    return seen_targets, extra_targets


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare the BARCODE rerun inputs for GoT-Multi by reusing processed_input from the MUT run."
    )
    parser.add_argument("--mutation-params", required=True, help="Path to the MUT/gene-sequence params file")
    parser.add_argument("--barcode-params", required=True, help="Path to the BARCODE/probe-barcode params file")
    parser.add_argument("--processed-input-dir", required=True, help="processed_input directory copied into the BARCODE output directory")
    parser.add_argument("--out-params", required=True, help="Path to write the generated BARCODE params file")
    args = parser.parse_args()

    mutation_params = Path(args.mutation_params).resolve()
    barcode_params = Path(args.barcode_params).resolve()
    processed_input_dir = Path(args.processed_input_dir).resolve()
    out_params = Path(args.out_params).resolve()

    for path in (mutation_params, barcode_params):
        if not path.is_file():
            raise ParamsError(f"Params file not found: {path}")
    if not processed_input_dir.is_dir():
        raise ParamsError(f"processed_input directory not found: {processed_input_dir}")

    mut_order, mut_sections = parse_params_file(mutation_params)
    _, bar_sections = parse_params_file(barcode_params)

    if "config-files" not in mut_sections or "config-files" not in bar_sections:
        raise ParamsError("Both params files must contain a [config-files] section")
    if "ironthrone-args" not in mut_sections or "ironthrone-args" not in bar_sections:
        raise ParamsError("Both params files must contain an [ironthrone-args] section")

    mut_args = parse_ironthrone_args(mut_sections["ironthrone-args"], mutation_params)
    bar_args = parse_ironthrone_args(bar_sections["ironthrone-args"], barcode_params)

    barcode_suffix = bar_args.get("suffix", "").strip()
    if not barcode_suffix:
        raise ParamsError("The BARCODE params file must define suffix in [ironthrone-args]")

    primary_prefix = mut_args.get("prefix", "").strip()
    secondary_prefix = bar_args.get("prefix", "").strip()
    prefixes = []
    for prefix in (primary_prefix, secondary_prefix):
        if prefix and prefix not in prefixes:
            prefixes.append(prefix)

    barcode_configs = expand_config_entries(bar_sections["config-files"], barcode_params)
    barcode_config_by_target: dict[str, Path] = {}
    for config_path in barcode_configs:
        target = extract_target_name(config_path, prefixes, barcode_suffix)
        if target in barcode_config_by_target:
            raise ParamsError(f"Duplicate BARCODE config target {target!r}: {config_path}")
        barcode_config_by_target[target] = config_path

    generated_sections = OrderedDict((section, list(lines)) for section, lines in mut_sections.items())
    generated_sections["config-files"] = [str(path) for path in barcode_configs]

    generated_args = OrderedDict(mut_args)
    generated_args["suffix"] = barcode_suffix
    generated_sections["ironthrone-args"] = [f"{key},{value}" for key, value in generated_args.items()]
    write_params_file(mut_order, generated_sections, out_params)

    processed_params_csv = find_processed_params_csv(processed_input_dir)
    seen_targets, _ = rewrite_processed_params_csv(processed_params_csv, barcode_config_by_target)

    print(f"Generated BARCODE params: {out_params}")
    print(f"Rewrote processed_input params CSV: {processed_params_csv}")
    print(f"Targets prepared for BARCODE rerun ({len(seen_targets)}): {', '.join(seen_targets)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ParamsError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
