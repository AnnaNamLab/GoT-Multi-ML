#!/usr/bin/env python3
"""
Process IronThrone Output (Paired geneseq, probebc) and produce ironthrone_out_pp.csv

This script processes IronThrone-classic output for paired mutation and barcode targets, merges with GEX data, and produces a unified CSV for IronThrone-ML.

Arguments:
    --mutdir: Directory containing mutation (geneseq) IronThrone results
    --bardir: Directory containing barcode (probebc) IronThrone results
    --cell-group: Column name in GEX data for cell group
    --negative-control-cell-group: Comma-separated list of negative control cell groups
    --gex: Path to GEX .h5ad file
    --target-info: Path to target info CSV
    --outdir: Output directory
    --additional-features: Comma-separated list of additional features to add from target info (optional)

Output:
    ironthrone_out_pp.csv in the output directory
"""

import os
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ironthrone_ml.pp import paired_pp_ironthrone_result
from ironthrone_ml import createFolder


def get_matching_dirs(dir_path, pattern="*.summTable.concat.umi_collapsed.txt"):
    import glob

    matching_dirs = []
    for item in os.listdir(dir_path):
        if item != "processed_input":
            full_path = os.path.join(dir_path, item)
            if os.path.isdir(full_path):
                if glob.glob(os.path.join(full_path, pattern)):
                    matching_dirs.append(item)
    return matching_dirs


def main():
    parser = argparse.ArgumentParser(
        description="Process IronThrone output and produce ironthrone_out_pp.csv"
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
    parser.add_argument(
        "--cell-group", type=str, required=True, help="Cell group column in GEX data"
    )
    parser.add_argument(
        "--negative-control-cell-group",
        type=str,
        default="",
        help="Comma-separated list of negative control cell groups",
    )
    parser.add_argument("--gex", type=str, required=True, help="Path to GEX .h5ad file")
    parser.add_argument(
        "--target-info", type=str, required=True, help="Path to target info CSV"
    )
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--additional-features",
        type=str,
        default="",
        help="Comma-separated list of additional features (optional)",
    )
    args = parser.parse_args()

    createFolder(args.outdir)
    negative_control_cell_group = [
        x for x in args.negative_control_cell_group.split(",") if x
    ]
    additional_features = [x for x in args.additional_features.split(",") if x]

    gex = sc.read_h5ad(args.gex)
    gex.obs["is_non_mut"] = gex.obs[args.cell_group].apply(
        lambda x: True if x in negative_control_cell_group else False
    )
    gex_df = gex.obs[
        ["BC", "experiment", "sample", "is_non_mut", args.cell_group]
    ].copy()

    mut_targets = get_matching_dirs(args.mutdir)
    bar_targets = get_matching_dirs(args.bardir)
    targets_of_interest = list(set(bar_targets) & set(mut_targets))
    print(f"Targets of interest: {targets_of_interest}")

    probe_info = pd.read_csv(args.target_info)
    gen_df = []
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
        target_gen_df = pd.merge(gex_df, gen, on="BC", how="inner")
        assert target_gen_df["BC"].nunique() == target_gen_df.shape[0]
        target_gen_df.set_index("BC", inplace=True)
        target_gen_df["expected_sample"] = probe_info.query(f'target == "{target}"')[
            "sample"
        ].values[0]
        gen_df.append(target_gen_df)
    if not gen_df:
        print("No valid targets found.")
        return
    gen_df = pd.concat(gen_df).reset_index()
    targets = gen_df.sort_values("expected_sample")["target"].unique()
    print("#Targets:", len(targets))
    gen_df["in_expected_sample"] = (
        gen_df["sample"] == gen_df["expected_sample"]
    ).astype(int)
    gen_df["Y"] = 1
    gen_df.loc[
        (gen_df["in_expected_sample"] == 1)
        & (gen_df["genotype_merged"] == "MUT")
        & (gen_df["is_non_mut"] == False),
        "Y",
    ] = 2
    gen_df.loc[
        (gen_df["in_expected_sample"] == 0) & (gen_df["genotype_merged"] == "MUT"), "Y"
    ] = 0
    gen_df.loc[
        (gen_df["genotype_merged"] == "MUT") & (gen_df["is_non_mut"] == True), "Y"
    ] = 0
    gen_df.loc[(gen_df["genotype_merged"] != "MUT"), "Y"] = 1
    int2label = {0: "FalsePositive", 1: "WT", 2: "MUT"}
    gen_df["Y_num"] = gen_df["Y"].copy()
    gen_df["Y"] = gen_df["Y"].map(int2label)
    if "is_non_mut" in gen_df.columns:
        gen_df["expected"] = (
            (gen_df["sample"] == gen_df["expected_sample"]) & (~gen_df["is_non_mut"])
        ).astype(int)
    else:
        gen_df["expected"] = (gen_df["sample"] == gen_df["expected_sample"]).astype(int)
    if additional_features:
        for feature in additional_features:
            if feature not in probe_info.columns:
                print(f"Feature {feature} not found in probe_info.")
                continue
            target2feature = (
                probe_info[["target", feature]]
                .drop_duplicates()
                .set_index("target")
                .to_dict()[feature]
            )
            if probe_info["target"].nunique() != len(target2feature):
                print(
                    f"Feature {feature} has duplicate values for targets. Skipping..."
                )
                continue
            gen_df[feature] = gen_df["target"].map(target2feature)
    unified_columns = [
        "BC",
        "target",
        "sample",
        "expected_sample",
        "expected",
        "experiment",
        "is_non_mut",
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
        "wt_calls_geneseq",
        "mut_calls_geneseq",
        "amb_calls_geneseq",
        "wt_calls_probebc",
        "mut_calls_probebc",
        "amb_calls_probebc",
        "genotype_merged",
        "Y",
        "Y_num",
    ]
    unified_columns += additional_features
    assert gen_df[unified_columns].isna().sum().sum() == 0
    gen_df = gen_df[unified_columns]
    gen_df.to_csv(os.path.join(args.outdir, "ironthrone_out_pp.csv"), index=False)
    print(f"Output written to {os.path.join(args.outdir, 'ironthrone_out_pp.csv')}")


if __name__ == "__main__":
    main()
