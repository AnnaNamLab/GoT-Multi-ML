#!/usr/bin/env python3
"""Prepare ML input features from the paired GoT-Multi genotype profile.

This script adds GEX metadata and target-level annotations to the paired
ironthrone_out_pp.csv produced by the non-ML GoT-Multi wrapper, and creates the
mandatory label columns used for optional ML denoising.
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils import createFolder

BASE_REQUIRED_GEX_COLUMNS = ["BC", "sample"]
RESERVED_TARGET_INFO_COLUMNS = {"target", "sample", "expected_sample", "experiment"}
LEGACY_DERIVED_COLUMNS = {
    "sample",
    "experiment",
    "cell_group",
    "is_non_mut",
    "expected_sample",
    "in_expected_sample",
    "expected",
    "Y",
    "Y_num",
}


TARGET_INFO_NON_FEATURE_COLUMNS = ["sample", "expected_sample", "experiment"]


def log_dropped_target_info_feature_columns(target_info: pd.DataFrame) -> None:
    present = [c for c in TARGET_INFO_NON_FEATURE_COLUMNS if c in target_info.columns]
    if not present:
        return

    print(
        "Dropping target_info columns from ML features: "
        f"{present}"
    )
    if any(c in present for c in ["sample", "expected_sample"]):
        print(
            "target_info sample/expected_sample columns are still used to define expected samples for labels."
        )


def parse_csv_list(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def parse_expected_sample_value(value) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(",") if x.strip()]


def normalize_bool_series(series: pd.Series, column_name: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(int).astype(bool)

    normalized = (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"true": True, "false": False, "1": True, "0": False})
    )
    if normalized.isna().any():
        bad_values = sorted(series[normalized.isna()].astype(str).unique().tolist())
        raise ValueError(
            f"Column {column_name} could not be interpreted as boolean. "
            f"Example unsupported values: {bad_values[:5]}"
        )
    return normalized


def build_gex_df(
    gex_obs: pd.DataFrame,
    cell_group_column: str,
    additional_gex_features: list[str],
    negative_control_cell_group: list[str],
) -> pd.DataFrame:
    required_gex_columns = BASE_REQUIRED_GEX_COLUMNS + [cell_group_column]
    missing_required = [c for c in required_gex_columns if c not in gex_obs.columns]
    if missing_required:
        raise ValueError(
            f"GEX object is missing required .obs columns: {missing_required}. "
            f"Expected at least {required_gex_columns}."
        )

    missing_additional = [c for c in additional_gex_features if c not in gex_obs.columns]
    if missing_additional:
        raise ValueError(
            f"Requested --additional-features were not found in GEX .obs: {missing_additional}"
        )

    selected_columns = BASE_REQUIRED_GEX_COLUMNS.copy()
    if "experiment" in gex_obs.columns:
        selected_columns.append("experiment")
    selected_columns.extend(
        [c for c in additional_gex_features if c not in selected_columns and c != "is_non_mut"]
    )

    gex_df = gex_obs[selected_columns].copy()
    gex_df["cell_group"] = gex_obs[cell_group_column].values
    if "experiment" not in gex_df.columns:
        print("GEX .obs does not contain 'experiment'. Filling with NA.")
        gex_df["experiment"] = pd.NA

    if "is_non_mut" in gex_obs.columns:
        if negative_control_cell_group:
            print(
                "GEX .obs already contains 'is_non_mut'; ignoring --negative-control-cell-group."
            )
        gex_df["is_non_mut"] = normalize_bool_series(gex_obs["is_non_mut"], "is_non_mut")
    else:
        if not negative_control_cell_group:
            raise ValueError(
                "GEX .obs is missing 'is_non_mut'. Provide --negative-control-cell-group "
                "to derive it from the configured cell-group column."
            )
        gex_df["is_non_mut"] = gex_df["cell_group"].isin(negative_control_cell_group)

    return gex_df


def aggregate_target_info(target_info: pd.DataFrame):
    if "target" not in target_info.columns:
        raise ValueError("target_info must contain a 'target' column.")

    if "expected_sample" in target_info.columns:
        expected_sample_source = "expected_sample"
    elif "sample" in target_info.columns:
        expected_sample_source = "sample"
    else:
        raise ValueError(
            "target_info must contain either 'expected_sample' or 'sample' so expected_sample can be derived."
        )

    target_info = target_info.drop_duplicates().copy()
    feature_columns = [
        c for c in target_info.columns if c not in RESERVED_TARGET_INFO_COLUMNS
    ]

    expected_rows = []
    feature_rows = []
    expected_sample_map = {}

    for target, group in target_info.groupby("target", sort=True):
        expected_samples = sorted(
            {
                expected_sample
                for value in group[expected_sample_source].tolist()
                for expected_sample in parse_expected_sample_value(value)
            }
        )
        if not expected_samples:
            raise ValueError(f"Target {target} does not have any expected sample values.")

        expected_rows.append(
            {"target": target, "expected_sample": ", ".join(expected_samples)}
        )
        expected_sample_map[target] = set(expected_samples)

        feature_row = {"target": target}
        for column in feature_columns:
            non_null_values = group[column].dropna().tolist()
            unique_values = pd.unique(non_null_values)
            if len(unique_values) > 1:
                preview = [str(x) for x in unique_values[:5]]
                raise ValueError(
                    f"target_info column {column!r} has multiple values for target {target!r}: {preview}. "
                    "Target-level ML features must be consistent across rows for the same target."
                )
            feature_row[column] = unique_values[0] if len(unique_values) == 1 else pd.NA
        feature_rows.append(feature_row)

    expected_df = pd.DataFrame(expected_rows)
    feature_df = pd.DataFrame(feature_rows)
    return expected_df, feature_df, expected_sample_map


def rename_target_feature_columns(
    feature_df: pd.DataFrame, existing_columns: set[str]
) -> tuple[pd.DataFrame, dict[str, str]]:
    rename_map = {}
    used_columns = set(existing_columns)
    for column in feature_df.columns:
        if column == "target":
            continue
        candidate = column
        if candidate in used_columns:
            candidate = f"target_info_{column}"
        while candidate in used_columns:
            candidate = f"target_info_{candidate}"
        rename_map[column] = candidate
        used_columns.add(candidate)

    return feature_df.rename(columns=rename_map), rename_map


def drop_legacy_derived_columns(
    pp_df: pd.DataFrame,
    additional_gex_features: list[str],
    target_info: pd.DataFrame,
) -> pd.DataFrame:
    candidate_columns = set(LEGACY_DERIVED_COLUMNS)
    candidate_columns.update(additional_gex_features)
    candidate_columns.update(c for c in target_info.columns if c != "target")
    candidate_columns.update(
        f"target_info_{c}"
        for c in target_info.columns
        if c not in RESERVED_TARGET_INFO_COLUMNS
    )

    drop_columns = [c for c in pp_df.columns if c in candidate_columns]
    if drop_columns:
        print(
            "Dropping legacy derived columns from input ironthrone_out_pp.csv: "
            f"{drop_columns}"
        )
        pp_df = pp_df.drop(columns=drop_columns)
    return pp_df


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare ML denoising input from paired GoT-Multi genotype output"
    )
    parser.add_argument(
        "--ironthrone-pp", required=True, help="Path to paired ironthrone_out_pp.csv"
    )
    parser.add_argument("--gex", required=True, help="Path to GEX .h5ad file")
    parser.add_argument(
        "--target-info", required=True, help="Path to target info CSV"
    )
    parser.add_argument(
        "--cell-group",
        default="cell_group",
        help="Column in GEX .obs to use as the cell-group annotation (default: cell_group)",
    )
    parser.add_argument(
        "--negative-control-cell-group",
        default="",
        help="Comma-separated negative-control cell groups. Used only if GEX .obs lacks is_non_mut.",
    )
    parser.add_argument(
        "--additional-features",
        default="",
        help="Comma-separated additional GEX .obs columns to include as ML features",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    createFolder(args.outdir)
    pp_df = pd.read_csv(args.ironthrone_pp)
    target_info = pd.read_csv(args.target_info)
    negative_control_cell_group = parse_csv_list(args.negative_control_cell_group)
    additional_gex_features = parse_csv_list(args.additional_features)

    pp_df = drop_legacy_derived_columns(pp_df, additional_gex_features, target_info)
    log_dropped_target_info_feature_columns(target_info)

    gex = sc.read_h5ad(args.gex)
    gex_df = build_gex_df(
        gex.obs,
        args.cell_group,
        additional_gex_features,
        negative_control_cell_group,
    )

    gen_barcodes = set(pp_df["BC"].dropna().astype(str))
    gex_barcodes = set(gex_df["BC"].dropna().astype(str))
    shared_barcodes = gen_barcodes & gex_barcodes
    print(
        "Intersecting GEN and GEX by BC: "
        f"GEN barcodes={len(gen_barcodes)}, "
        f"GEX barcodes={len(gex_barcodes)}, "
        f"shared={len(shared_barcodes)}"
    )

    merged_df = pd.merge(pp_df, gex_df, on="BC", how="inner")
    assert merged_df["BC"].astype(str).isin(gex_barcodes).all()
    print(
        "Filtered GEN rows to cells present in GEX: "
        f"{len(pp_df)} -> {len(merged_df)} rows"
    )

    expected_df, target_feature_df, expected_sample_map = aggregate_target_info(target_info)
    merged_df = pd.merge(merged_df, expected_df, on="target", how="left")
    if merged_df["expected_sample"].isna().any():
        missing_targets = sorted(
            merged_df.loc[merged_df["expected_sample"].isna(), "target"].unique()
        )
        raise ValueError(f"Missing target_info mapping for targets: {missing_targets}")

    target_feature_df, rename_map = rename_target_feature_columns(
        target_feature_df, set(merged_df.columns)
    )
    merged_df = pd.merge(merged_df, target_feature_df, on="target", how="left")
    merged_df["in_expected_sample"] = merged_df.apply(
        lambda row: str(row["sample"]).strip() in expected_sample_map.get(row["target"], set()),
        axis=1,
    ).astype(int)
    merged_df["expected"] = (
        merged_df["in_expected_sample"].astype(bool) & (~merged_df["is_non_mut"])
    ).astype(int)

    merged_df["Y_num"] = 1
    merged_df.loc[
        (merged_df["in_expected_sample"] == 1)
        & (merged_df["genotype_merged"] == "MUT")
        & (~merged_df["is_non_mut"]),
        "Y_num",
    ] = 2
    merged_df.loc[
        (merged_df["in_expected_sample"] == 0)
        & (merged_df["genotype_merged"] == "MUT"),
        "Y_num",
    ] = 0
    merged_df.loc[
        (merged_df["genotype_merged"] == "MUT") & (merged_df["is_non_mut"]),
        "Y_num",
    ] = 0
    merged_df.loc[merged_df["genotype_merged"] != "MUT", "Y_num"] = 1
    merged_df["Y"] = merged_df["Y_num"].map({0: "FalsePositive", 1: "WT", 2: "MUT"})

    if rename_map:
        print(f"Included target_info feature columns: {list(rename_map.values())}")
    if additional_gex_features:
        print(f"Included additional GEX feature columns: {additional_gex_features}")

    output_path = os.path.join(args.outdir, "ironthrone_out_pp_ml_input.csv")
    merged_df.to_csv(output_path, index=False)
    print(f"ML input written to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
