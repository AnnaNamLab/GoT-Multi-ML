#!/usr/bin/env python3
"""Prepare the MUT first pass for GoT-Multi.

This helper enforces the GoT-Multi rule that expected sequences must come from
mutation / gene-sequence inputs only. If the MUT params file already specifies
an expected-sequence CSV, that CSV is reused. Otherwise a new CSV is generated
from the MUT configs and the generated params file is updated to point to it.
"""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from pathlib import Path

from build_expected_sequence_csv import write_expected_sequence_csv
from prepare_got_multi_barcode_run import (
    ParamsError,
    expand_config_entries,
    parse_ironthrone_args,
    parse_params_file,
    write_params_file,
)


def resolve_optional_path(raw_path: str, source: Path) -> Path:
    candidate = Path(raw_path).expanduser()
    if not candidate.is_absolute():
        candidate = source.parent / candidate
    return candidate.resolve()


def get_expected_sequence_csv_from_params(
    sections: OrderedDict[str, list[str]],
    source: Path,
) -> Path | None:
    lines = [line.strip() for line in sections.get("expected-sequence-csv", []) if line.strip()]
    if not lines:
        return None
    if len(lines) > 1:
        raise ParamsError(
            f"[expected-sequence-csv] must contain at most one path in {source}, found {len(lines)} entries"
        )

    expected_sequence_csv = resolve_optional_path(lines[0], source)
    if not expected_sequence_csv.is_file():
        raise ParamsError(
            f"Expected-sequence CSV specified in {source} was not found: {expected_sequence_csv}"
        )
    return expected_sequence_csv


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the MUT first-pass inputs for GoT-Multi. Reuse the expected-sequence CSV "
            "from the MUT params file when present; otherwise generate it from the MUT configs."
        )
    )
    parser.add_argument("--mutation-params", required=True, help="Path to the MUT/gene-sequence params file")
    parser.add_argument("--out-params", required=True, help="Path to write the generated MUT params file")
    parser.add_argument(
        "--out-expected-sequence-csv",
        required=True,
        help="Path to write the generated expected-sequence CSV if one is not already specified",
    )
    args = parser.parse_args()

    mutation_params = Path(args.mutation_params).resolve()
    out_params = Path(args.out_params).resolve()
    out_expected_sequence_csv = Path(args.out_expected_sequence_csv).resolve()

    if not mutation_params.is_file():
        raise ParamsError(f"Params file not found: {mutation_params}")

    order, sections = parse_params_file(mutation_params)
    if "config-files" not in sections:
        raise ParamsError(f"[config-files] section is missing in {mutation_params}")
    if "ironthrone-args" not in sections:
        raise ParamsError(f"[ironthrone-args] section is missing in {mutation_params}")

    args_map = parse_ironthrone_args(sections["ironthrone-args"], mutation_params)
    prefix = args_map.get("prefix", "").strip()
    suffix = args_map.get("suffix", ".config").strip() or ".config"
    mutation_configs = expand_config_entries(sections["config-files"], mutation_params)

    expected_sequence_csv = get_expected_sequence_csv_from_params(sections, mutation_params)
    expected_sequence_mode = "provided"
    if expected_sequence_csv is None:
        write_expected_sequence_csv(
            mutation_configs,
            output=out_expected_sequence_csv,
            prefix=prefix,
            suffix=suffix,
        )
        expected_sequence_csv = out_expected_sequence_csv
        expected_sequence_mode = "generated"

    generated_sections = OrderedDict((section, list(lines)) for section, lines in sections.items())
    had_expected_sequence_section = "expected-sequence-csv" in generated_sections
    generated_sections["expected-sequence-csv"] = [str(expected_sequence_csv)]

    if not had_expected_sequence_section:
        if "config-files" in order:
            insert_at = order.index("config-files") + 1
        else:
            insert_at = len(order)
        order = order[:insert_at] + ["expected-sequence-csv"] + order[insert_at:]

    write_params_file(order, generated_sections, out_params)
    print(f"Generated MUT params: {out_params}")
    if expected_sequence_mode == "generated":
        print(f"Generated expected-sequence CSV from MUT configs: {expected_sequence_csv}")
    else:
        print(f"Reusing expected-sequence CSV from MUT params: {expected_sequence_csv}")
    print(f"MUT targets prepared ({len(mutation_configs)}): {', '.join(path.name for path in mutation_configs)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ParamsError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
