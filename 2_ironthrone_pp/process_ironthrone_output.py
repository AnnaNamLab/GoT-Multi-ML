#!/usr/bin/env python3
"""Process paired IronThrone-Multi outputs and produce a merged genotype profile.

This step merges the GeneSeq and ProbeBC IronThrone results target-by-target,
creates genotype_geneseq, genotype_probebc, and genotype_merged, and writes the
combined GoT-Multi paired genotype profile to ironthrone_out_pp.csv.

This stage is intentionally GEX-free. GEX metadata and ML labels are added later
by the optional ML preparation / denoising stage.
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils.pp import paired_pp_ironthrone_result
from utils import createFolder


READ_STAT_COLUMNS = [
    "wt_reads_per_umi_avg_geneseq",
    "wt_reads_per_umi_med_geneseq",
    "wt_reads_per_umi_std_geneseq",
    "wt_reads_per_umi_total_geneseq",
    "wt_reads_per_umi_count_geneseq",
    "wt_reads_per_umi_avg_probebc",
    "wt_reads_per_umi_med_probebc",
    "wt_reads_per_umi_std_probebc",
    "wt_reads_per_umi_total_probebc",
    "wt_reads_per_umi_count_probebc",
    "mut_reads_per_umi_avg_geneseq",
    "mut_reads_per_umi_med_geneseq",
    "mut_reads_per_umi_std_geneseq",
    "mut_reads_per_umi_total_geneseq",
    "mut_reads_per_umi_count_geneseq",
    "mut_reads_per_umi_avg_probebc",
    "mut_reads_per_umi_med_probebc",
    "mut_reads_per_umi_std_probebc",
    "mut_reads_per_umi_total_probebc",
    "mut_reads_per_umi_count_probebc",
    "amb_reads_per_umi_avg_geneseq",
    "amb_reads_per_umi_med_geneseq",
    "amb_reads_per_umi_std_geneseq",
    "amb_reads_per_umi_total_geneseq",
    "amb_reads_per_umi_count_geneseq",
    "amb_reads_per_umi_avg_probebc",
    "amb_reads_per_umi_med_probebc",
    "amb_reads_per_umi_std_probebc",
    "amb_reads_per_umi_total_probebc",
    "amb_reads_per_umi_count_probebc",
]

CALL_COLUMNS = [
    "wt_calls_geneseq",
    "mut_calls_geneseq",
    "amb_calls_geneseq",
    "wt_calls_probebc",
    "mut_calls_probebc",
    "amb_calls_probebc",
]

GENOTYPE_COLUMNS = [
    "genotype_geneseq",
    "genotype_probebc",
    "genotype_merged",
]

UNIFIED_COLUMNS = ["BC", "target"] + READ_STAT_COLUMNS + CALL_COLUMNS + GENOTYPE_COLUMNS


def get_matching_dirs(dir_path: str, pattern: str = "*.summTable.concat.umi_collapsed.txt") -> list[str]:
    matching_dirs: list[str] = []
    for item in os.listdir(dir_path):
        if item == "processed_input":
            continue
        full_path = os.path.join(dir_path, item)
        if os.path.isdir(full_path) and glob.glob(os.path.join(full_path, pattern)):
            matching_dirs.append(item)
    return matching_dirs


def write_stage1_genotype_count_plot(gen_df: pd.DataFrame, outdir: str) -> str:
    genotype_order = ["WT", "MUT", "Ambiguous", "Unprofiled"]
    palette = {
        "WT": "#8ecae6",
        "MUT": "#023e8a",
        "Ambiguous": "#bdbdbd",
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
        counts = target_df["genotype_merged"].value_counts()
        x_values = genotype_order
        y_values = [int(counts.get(g, 0)) for g in x_values]
        colors = [palette[g] for g in x_values]

        bars = ax.bar(x_values, y_values, color=colors, edgecolor="0.2", linewidth=0.8)
        ax.set_title(target, fontsize=9, weight="bold")
        ax.set_xlabel("genotype_merged", fontsize=8)
        ax.set_ylabel("cell count", fontsize=8)
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
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
    output_path = os.path.join(outdir, "initial_genotype_counts_by_target.pdf")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Stage-1 genotype count plot written to {output_path}")
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge paired GeneSeq and ProbeBC IronThrone outputs into ironthrone_out_pp.csv"
    )
    parser.add_argument(
        "--mutdir",
        type=str,
        required=True,
        help="Directory with mutation (geneseq) IronThrone results",
    )
    parser.add_argument(
        "--bardir",
        type=str,
        required=True,
        help="Directory with barcode (probebc) IronThrone results",
    )
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    args = parser.parse_args()

    createFolder(args.outdir)

    mut_targets = set(get_matching_dirs(args.mutdir))
    bar_targets = set(get_matching_dirs(args.bardir))
    targets_of_interest = sorted(mut_targets & bar_targets)
    print(f"Targets of interest: {targets_of_interest}")

    merged_frames: list[pd.DataFrame] = []
    for target in targets_of_interest:
        gen = paired_pp_ironthrone_result(
            geneseq_dir=args.mutdir,
            probebc_dir=args.bardir,
            target=target,
            calculate_reads_stats=True,
        )
        if gen is None:
            print(f"Skipping {target} due to missing data.")
            continue

        missing_columns = [column for column in UNIFIED_COLUMNS if column not in gen.columns]
        if missing_columns:
            raise ValueError(f"Target {target} is missing required columns: {missing_columns}")

        target_gen_df = gen[UNIFIED_COLUMNS].copy()
        if target_gen_df["BC"].duplicated().any():
            duplicated_bc = target_gen_df.loc[target_gen_df["BC"].duplicated(), "BC"].unique()[:5]
            raise ValueError(
                f"Target {target} does not have unique BC rows after paired post-processing. "
                f"Example duplicated BCs: {duplicated_bc.tolist()}"
            )
        merged_frames.append(target_gen_df)

    if not merged_frames:
        print("No valid targets found.")
        return 0

    gen_df = pd.concat(merged_frames, ignore_index=True)
    gen_df = gen_df.sort_values(["target", "BC"]).reset_index(drop=True)

    output_path = os.path.join(args.outdir, "ironthrone_out_pp.csv")
    gen_df.to_csv(output_path, index=False)
    print(f"Output written to {output_path}")

    write_stage1_genotype_count_plot(gen_df, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
