#!/usr/bin/env python3
"""
Post-process single IronThrone-Multi outputs into per-target and combined CSV tables.

This script handles one IronThrone-Multi run directory at a time.
It does not pair mutation/barcode runs or merge with GEX metadata. Instead it:

1. Reads each `<TARGET>.summTable.concat.umi_collapsed.txt`
2. Renames columns to the standardized single-run output naming scheme
3. Computes per-cell read-count summary statistics
4. Assigns a single-run genotype label
5. Writes:
   - `<TARGET>/<TARGET>.ironthrone_out_pp.csv` for each target
   - `<RESULTDIR>/ironthrone_out_pp.csv` combining all targets
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RENAME_MAP = {
    "num.WT.in.dups": "wt_reads_per_umi",
    "num.MUT.in.dups": "mut_reads_per_umi",
    "num.amb.in.dups": "amb_reads_per_umi",
    "call.in.dups": "calls_per_umi",
    "WT.calls": "wt_calls",
    "MUT.calls": "mut_calls",
    "amb.calls": "amb_calls",
}

READ_COUNT_COLUMNS = [
    "wt_reads_per_umi",
    "mut_reads_per_umi",
    "amb_reads_per_umi",
]

CORE_COLUMNS = [
    "BC",
    "target",
    "UMI",
    "wt_reads_per_umi",
    "mut_reads_per_umi",
    "amb_reads_per_umi",
    "calls_per_umi",
    "wt_calls",
    "mut_calls",
    "amb_calls",
    "genotype",
]


def split_numeric_values(value: object) -> List[float]:
    if pd.isna(value):
        return [0.0]

    values: List[float] = []
    for item in str(value).split(";"):
        item = item.strip()
        if item in {"", "nan", "NaN", "None"}:
            continue
        try:
            values.append(float(item))
        except ValueError:
            continue

    if not values:
        return [0.0]
    return values


def stat_columns(reads_col: str) -> List[str]:
    return [
        f"{reads_col}_avg",
        f"{reads_col}_med",
        f"{reads_col}_std",
        f"{reads_col}_total",
        f"{reads_col}_count",
    ]


ALL_OUTPUT_COLUMNS = list(CORE_COLUMNS)
for _reads_col in READ_COUNT_COLUMNS:
    ALL_OUTPUT_COLUMNS.extend(stat_columns(_reads_col))


def ordered_columns(df: pd.DataFrame) -> List[str]:
    return [col for col in ALL_OUTPUT_COLUMNS if col in df.columns]


def target_result_path(resultdir: Path, target: str) -> Path:
    return resultdir / target / f"{target}.summTable.concat.umi_collapsed.txt"


def discover_targets(resultdir: Path) -> List[str]:
    targets: List[str] = []
    for item in sorted(resultdir.iterdir()):
        if not item.is_dir() or item.name == "processed_input":
            continue
        if target_result_path(resultdir, item.name).exists():
            targets.append(item.name)
    return targets


def compute_read_stats(df: pd.DataFrame) -> None:
    for reads_col in READ_COUNT_COLUMNS:
        split_values = df[reads_col].apply(split_numeric_values)
        df[f"{reads_col}_avg"] = split_values.apply(np.mean)
        df[f"{reads_col}_med"] = split_values.apply(np.median)
        df[f"{reads_col}_std"] = split_values.apply(np.std)
        df[f"{reads_col}_total"] = split_values.apply(np.sum)
        df[f"{reads_col}_count"] = split_values.apply(len).astype(int)


def pp_ironthrone_result(
    resultdir: Path,
    target: str,
    calculate_reads_stats: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    result_file = target_result_path(resultdir, target)
    df = (
        pd.read_csv(result_file, sep="\t")
        .sort_values(["BC", "UMI"], kind="stable")
        .rename(columns=RENAME_MAP)
    )

    for col in ["wt_calls", "mut_calls", "amb_calls"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    df["target"] = target
    df["genotype"] = "Unprofiled"
    df.loc[df["wt_calls"] > 0, "genotype"] = "WT"
    df.loc[df["mut_calls"] > 0, "genotype"] = "MUT"

    if calculate_reads_stats:
        compute_read_stats(df)

    df = df[ordered_columns(df)]

    if verbose:
        print(f"genotype for {target}:")
        print(df["genotype"].value_counts(dropna=False))
        print("")

    return df


def write_csv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, float_format="%.17g")


def write_single_run_genotype_count_plot(gen_df: pd.DataFrame, outdir: Path) -> Path:
    genotype_order = ["WT", "MUT", "Unprofiled"]
    palette = {
        "WT": "#8ecae6",
        "MUT": "#023e8a",
        "Unprofiled": "#e0e0e0",
    }

    targets = gen_df["target"].dropna().unique().tolist()
    n_targets = len(targets)
    n_cols = 2 if n_targets > 1 else 1
    n_rows = math.ceil(n_targets / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4.0, n_rows * 3.4))
    try:
        axes = axes.flatten()
    except AttributeError:
        axes = [axes]

    for i, target in enumerate(targets):
        ax = axes[i]
        target_df = gen_df[gen_df["target"] == target].copy()
        counts = target_df["genotype"].value_counts()
        x_values = genotype_order
        y_values = [int(counts.get(g, 0)) for g in x_values]
        colors = [palette[g] for g in x_values]

        bars = ax.bar(x_values, y_values, color=colors, edgecolor="0.2", linewidth=0.8)
        ax.set_title(target, fontsize=9, weight="bold")
        ax.set_xlabel("genotype", fontsize=8)
        ax.set_ylabel("cell count", fontsize=8)
        ax.tick_params(axis="x", labelrotation=25, labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ymax = max(y_values) if y_values else 0
        ax.set_ylim(0, ymax * 1.18 + 1)
        for bar, value in zip(bars, y_values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(ymax * 0.02, 0.3),
                str(value),
                ha="center",
                va="bottom",
                fontsize=8,
                weight="bold",
                color="0.2",
            )

    for ax in axes[n_targets:]:
        ax.axis("off")

    fig.tight_layout()
    output_path = outdir / "initial_genotype_counts_by_target.pdf"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Single-run genotype count plot written to {output_path}")
    return output_path


def combine_targets(resultdir: Path, targets: Iterable[str], verbose: bool = False) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for target in targets:
        df = pp_ironthrone_result(resultdir, target, calculate_reads_stats=True, verbose=verbose)
        write_csv(resultdir / target / f"{target}.ironthrone_out_pp.csv", df)
        frames.append(df)

    if not frames:
        empty = pd.DataFrame(columns=ALL_OUTPUT_COLUMNS)
        return empty

    combined = pd.concat(frames, ignore_index=True)
    return combined.sort_values(["target", "BC", "UMI"], kind="stable").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Post-process single IronThrone-Multi outputs into ironthrone_out_pp.csv."
    )
    parser.add_argument(
        "--resultdir",
        required=True,
        type=Path,
        help="IronThrone-Multi result directory containing target subdirectories.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Optional output directory for the combined ironthrone_out_pp.csv (default: resultdir).",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=None,
        help="Optional target name to process. Can be passed multiple times.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print genotype counts while processing targets.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    resultdir = args.resultdir.resolve()
    outdir = (args.outdir or resultdir).resolve()

    if not resultdir.is_dir():
        raise FileNotFoundError(f"Result directory not found: {resultdir}")

    targets = list(args.target) if args.target else discover_targets(resultdir)
    if not targets:
        combined = pd.DataFrame(columns=ALL_OUTPUT_COLUMNS)
        write_csv(outdir / "ironthrone_out_pp.csv", combined)
        print(f"No targets found. Wrote empty output to {outdir / 'ironthrone_out_pp.csv'}")
        return 0

    combined = combine_targets(resultdir, targets, verbose=args.verbose)
    write_csv(outdir / "ironthrone_out_pp.csv", combined)
    print(f"Output written to {outdir / 'ironthrone_out_pp.csv'}")
    write_single_run_genotype_count_plot(combined, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
