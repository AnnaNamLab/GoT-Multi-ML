import os
import numpy as np
import pandas as pd


def pp_ironthrone_result(resultdir, target, calculate_reads_stats=True, verbose=False):
    """Import IronThrone results for a target and preprocess it.
    Args:
        resultdir: ironthrone result directory
        target: genotype target name (subdirectory name within resultdir)
        verbose (bool, optional): Prints out genotype counts. Defaults to False.
    Steps:
        1. Import IronThrone result file
        2. Rename some columns to be more informative and unified (to me)
        3. Add a 'genotype' column based on the calls
        4. if `calculate_reads_stats` is True, calculate some statistics for the read counts
    """
    ironthrone_result = (
        pd.read_csv(
            os.path.join(
                resultdir, target, f"{target}.summTable.concat.umi_collapsed.txt"
            ),
            sep="\t",
        )
        .sort_values(["BC", "UMI"])
        .rename(
            columns={
                "num.WT.in.dups": "wt_reads_per_umi",  # number of wt reads per umi (separated by ;)
                "num.MUT.in.dups": "mut_reads_per_umi",  # number of mut reads per umi (separated by ;)
                "num.amb.in.dups": "amb_reads_per_umi",  # number of ambiguous reads per umi (separated by ;)
                "call.in.dups": "calls_per_umi",  # genotype call (based on reads) per umi (separated by ;)
                "WT.calls": "wt_calls",  # number of wt umi calls of the cell
                "MUT.calls": "mut_calls",  # number of mut umi calls of the cell
                "amb.calls": "amb_calls",  # number of ambiguous umi calls of the cell
            }
        )
    )

    # 'MUT' if mut_calls > 0, 'WT' elif wt_calls > 0, 'Unprofiled' otherwise
    ironthrone_result["genotype"] = "Unprofiled"
    ironthrone_result.loc[ironthrone_result["wt_calls"] > 0, "genotype"] = "WT"
    ironthrone_result.loc[ironthrone_result["mut_calls"] > 0, "genotype"] = "MUT"

    # Calculate some statistics for the read counts
    if calculate_reads_stats is True:
        # * Calculate Reads statistics
        read_count_columns = [
            "wt_reads_per_umi",
            "mut_reads_per_umi",
            "amb_reads_per_umi",
        ]

        # Loop through each column and reads_col to create the new columns
        for reads_col in read_count_columns:
            ironthrone_result[f"{reads_col}_avg"] = (
                ironthrone_result[reads_col]
                .fillna(0)
                .astype(str)
                .str.split(";")
                .apply(lambda x: np.mean([float(i) for i in x if i != "nan"]))
            )
            ironthrone_result[f"{reads_col}_med"] = (
                ironthrone_result[reads_col]
                .fillna(0)
                .astype(str)
                .str.split(";")
                .apply(lambda x: np.median([float(i) for i in x if i != "nan"]))
            )
            ironthrone_result[f"{reads_col}_std"] = (
                ironthrone_result[reads_col]
                .fillna(0)
                .astype(str)
                .str.split(";")
                .apply(lambda x: np.std([float(i) for i in x if i != "nan"]))
            )
            ironthrone_result[f"{reads_col}_total"] = (
                ironthrone_result[reads_col]
                .fillna(0)
                .astype(str)
                .str.split(";")
                .apply(lambda x: np.sum([float(i) for i in x if i != "nan"]))
            )
            ironthrone_result[f"{reads_col}_count"] = (
                ironthrone_result[reads_col]
                .fillna(0)
                .astype(str)
                .str.split(";")
                .apply(lambda x: len([float(i) for i in x if i != "nan"]))
            )

    if verbose:
        print(f"genotype for {target}:")
        print(ironthrone_result["genotype"].value_counts())
        print("")
    return ironthrone_result


def paired_pp_ironthrone_result(geneseq_dir, probebc_dir, target, **kwargs):
    """Import and preprocess IronThrone results for a target (for both geneseq and probebc). - Suitable for FRP project data
    Args:
        geneseq_dir: ironthrone result directory for geneseq (mutation)
        probebc_dir: ironthrone result directory for probebc (barcode)
        target: genotype target name (subdirectory name within result directory)
    Steps:
        1. Import and propcess ironthrone results for geneseq
            1-1) Import IronThrone result file
            1-2) Rename columns to unified convention (according to me)
            1-3) Add a 'genotype' column based on the WT/MUT/amb calls
            1-4) if `calculate_reads_stats` is True, calculate some statistics for the read counts
        2. Import and process ironthrone results for probebc
            2) Same as geneseq
        3. Combine the genotype of geneseq and probebc for each cell and define "genotype_merged"
        4. Filter out cells with: no read, no call in both geneseq and probebc
    """
    print(f"Processing {target}...")

    # * === Import & Genotype assignment for IronThrone result of target (Mutation / geneseq) === * #
    gen_mut = pp_ironthrone_result(geneseq_dir, target, **kwargs)

    # * === Import & Genotype assignment for IronThrone result of target (Barcode / probebc) === * #
    gen_bar = pp_ironthrone_result(probebc_dir, target, **kwargs)

    # * Merge the Mutation and Barcode results
    if (
        gen_bar["BC"].isna().all()
        or gen_bar["UMI"].isna().all()
        or gen_mut["BC"].isna().all()
        or gen_mut["UMI"].isna().all()
    ):
        print(f"{target} has no information on BC or UMI.")
        return None

    gen = pd.merge(
        gen_mut,
        gen_bar,
        on=["BC", "UMI"],
        suffixes=("_geneseq", "_probebc"),
        how="inner",
    )
    gen["target"] = target

    # * === Genotype Merged === #
    # - 'MUT' if both geneseq_genotype and probebc_genotype are 'MUT'
    # - 'WT' if both geneseq_genotype and probebc_genotype are 'WT'
    # - 'Unprofiled' if both geneseq_genotype and probebc_genotype are 'Unprofiled'
    # - 'Ambiguous' if geneseq_genotype and probebc_genotype are different
    gen["genotype_merged"] = ""
    gen.loc[
        (gen["genotype_geneseq"] == "Unprofiled")
        & (gen["genotype_probebc"] == "Unprofiled"),
        "genotype_merged",
    ] = "Unprofiled"
    gen.loc[
        (gen["genotype_geneseq"] != gen["genotype_probebc"]),
        "genotype_merged",
    ] = "Ambiguous"
    gen.loc[
        (gen["genotype_geneseq"] == "WT") & (gen["genotype_probebc"] == "WT"),
        "genotype_merged",
    ] = "WT"
    gen.loc[
        (gen["genotype_geneseq"] == "MUT") & (gen["genotype_probebc"] == "MUT"),
        "genotype_merged",
    ] = "MUT"
    assert (gen["genotype_merged"] == "").sum() == 0

    # * === Filter out cells not genotyped (with: no read, no call in both geneseq and probebc) === * #
    # Adjust genotyped_merged
    ## if row has all these columns as NaN, set 'genotype_merged' to NaN
    reads_cols = gen.filter(
        regex="reads_per_umi_(geneseq|probebc)$"
    ).columns  # number of reads per umi (separated by ;)
    umi_calls_cols = [
        "calls_per_umi_geneseq",
        "calls_per_umi_probebc",
    ]  # genotype call (based on reads) per umi (separated by ;)
    bc_calls_cols = gen.columns[
        gen.columns.str.contains("calls_probebc")
    ]  # number of wt/mut/amb umi calls of the cell (probebc)
    seq_calls_cols = gen.columns[
        gen.columns.str.contains("calls_geneseq")
    ]  # number of wt/mut/amb umi calls of the cell (geneseq)

    gen["genotype_merged"] = gen.apply(
        lambda x: np.nan
        if x[reads_cols].isna().all()
        and x[umi_calls_cols].isna().all()
        and x[bc_calls_cols].sum() == 0
        and x[seq_calls_cols].sum() == 0
        and x["genotype_merged"] not in ["WT", "MUT"]
        else x["genotype_merged"],
        axis=1,
    )

    # Fill in NaN values in the reads_per_umi columns
    gen[reads_cols] = gen[reads_cols].fillna("0")

    # ! Filter out cells with genotype_merged == NaN
    print(f"Before removing genotype_merged == NaN: {gen.shape}")
    gen.dropna(subset=["genotype_merged"], inplace=True)
    print(f"After removing genotype_merged == NaN: {gen.shape}")

    return gen
