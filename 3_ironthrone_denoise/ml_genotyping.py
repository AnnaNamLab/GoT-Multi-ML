#!/usr/bin/env python3
"""
ML-based genotyping denoising (IronThrone-ML) pipeline.

This script takes the output from IronThrone-PP (ironthrone_out_pp.csv),
trains and evaluates ML models for denoising genotyping calls, and outputs
refined predictions for each target.

Arguments:
    --input-features: Path to ironthrone_out_pp.csv
    --outdir: Output directory
    --max-allowed-fpr: Maximum allowed false positive rate (default: 0.03)
    --alpha: Weighting parameter for custom performance metric (default: 0.8)

Output:
    - Per-target ML results and logs
    - Aggregated and final prediction CSVs in <outdir>/aggregated/
"""

import os
import numpy as np
import pandas as pd
import argparse
import warnings
import sys
import mw_ml as ml
import mw_ml.classification as mcl
from matplotlib.gridspec import GridSpec
from sklearn.model_selection import train_test_split

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from ironthrone_ml import createFolder, logging
from ironthrone_ml.classification import custom_performance_metric
import ironthrone_ml.vz as mvz
import ironthrone_ml.custom as mc
from ironthrone_ml.custom import evaluate_prediction

random_seed = 1024
np.random.seed(random_seed)
warnings.filterwarnings("ignore")

columns_to_drop = [
    "BC",
    "target",
    "sample",
    "expected_sample",
    "expected",
    "experiment",
    "is_non_mut",
    "genotype_merged",
]
int2label = {0: "FalsePositive", 1: "WT", 2: "MUT"}


def main():
    parser = argparse.ArgumentParser(
        description="ML-based genotyping refinement for IronThrone-ML pipeline"
    )
    parser.add_argument(
        "--input-features",
        type=str,
        required=True,
        help="Path to ironthrone_out_pp.csv",
    )
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--max-allowed-fpr",
        type=float,
        default=0.03,
        help="Maximum allowed false positive rate",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.8, help="Alpha for custom performance metric"
    )
    parser.add_argument(
        "--models-to-run",
        type=str,
        default="logistic_regression,random_forest,knn,naive_bayes,xgboost,mlp,gradient_boosting,hist_gradient_boosting,adaboost",
        help="Comma-separated list of models to run (default: logistic_regression,random_forest,knn,naive_bayes,xgboost,mlp,gradient_boosting,hist_gradient_boosting,adaboost)",
    )
    args = parser.parse_args()

    createFolder(args.outdir)
    total_df = pd.read_csv(args.input_features)
    assert (
        total_df.query('genotype_merged == "MUT"')["Y_num"]
        .value_counts()
        .index.isin([0, 2])
        .all()
    )
    total_df.drop("Y", axis=1, inplace=True)
    total_df.rename(columns={"Y_num": "Y"}, inplace=True)

    targets = total_df["target"].unique()
    print(f"#Targets: {len(targets)}")
    skipped_targets = []
    models_to_run = [x.strip() for x in args.models_to_run.split(",") if x.strip()]
    for target in targets:
        target_df = total_df[total_df["target"] == target]
        non_target_df = total_df[total_df["target"] != target]
        if target_df["expected_sample"].unique()[0] not in target_df["sample"].unique():
            print(f"Skipping {target} due to lack of cells in expected sample")
            skipped_targets.append(target)
            continue
        if (target_df["Y"] == 2).sum() <= 1:
            print(f"Skipping {target} due to lack of MUT in expected sample")
            skipped_targets.append(target)
            continue
        resultdir = os.path.join(args.outdir, target)
        createFolder(resultdir)
        logging(f"# === Target: {target} === #", resultdir, "log.txt")
        logging(target_df["genotype_merged"].value_counts(), resultdir, "log.txt")

        # * ================================================================================ * #
        # * === First Run (Train: Half of target + Other targets & Test: The other half) === * #
        # * ================================================================================ * #
        # Duplicate rare classes if needed
        if (target_df["Y"] == 0).sum() == 1:
            target_df = pd.concat([target_df, target_df[target_df["Y"] == 0]])
        if (target_df["Y"] == 1).sum() == 1:
            target_df = pd.concat([target_df, target_df[target_df["Y"] == 1]])
        target_df_1, target_df_2 = train_test_split(
            target_df, test_size=0.5, stratify=target_df["Y"], random_state=random_seed
        )
        logging(
            f"First Run: Train: {target_df_1.shape}, Test: {target_df_2.shape}",
            resultdir,
            "log.txt",
        )

        # Train: Half of target + Other targets, Test: The other half
        train_df = pd.concat([target_df_1, non_target_df])
        test_df = target_df_2.copy()
        logging("Removing Unprofiled cells from training data", resultdir, "log.txt")
        train_df = train_df.query('genotype_merged != "Unprofiled"')
        logging("Train set Labels before balancing by label:", resultdir, "log.txt")
        logging(train_df["Y"].value_counts(), resultdir, "log.txt")

        # Balance by Label (WT, MUT, False Positive)
        train_df = ml.balance_by_label(
            train_df,
            "Y",
            balance_values=[0, 2],
            min_freq=100,
            random_seed=random_seed,
            visualize=False,
        )
        logging("Train set Labels after balancing by label:", resultdir, "log.txt")
        logging(train_df["Y"].value_counts(), resultdir, "log.txt")
        train_df.to_csv(os.path.join(resultdir, f"{target}_1_train.csv"), index=False)
        test_df.to_csv(os.path.join(resultdir, f"{target}_1_test.csv"), index=False)

        # Drop columns that are not features
        train_df = train_df.drop(columns=columns_to_drop)
        test_df = test_df.drop(columns=columns_to_drop)
        (
            test_res_1,
            fitted_models,
            _best_model_name,
            _y_pred_best,
            _y_pred_prob_best,
            _best_clf,
        ) = mcl.train_classification_model(
            train_df,
            test_df,
            resultdir,
            target_column="Y",
            y_value_of_interest=2,
            models_to_run=models_to_run,
            ensemble_hard=True,
            ensemble_soft=True,
            random_seed=random_seed,
        )

        # * ================================================================================= * #
        # * === Second Run (Train: Half of target + Other targets & Test: The other half) === * #
        # * ================================================================================= * #
        logging(
            f"Second Run: Train: {target_df_2.shape}, Test: {target_df_1.shape}",
            resultdir,
            "log.txt",
        )

        # Train: Half of target + Other targets, Test: The other half
        train_df = pd.concat([target_df_2, non_target_df])
        test_df = target_df_1.copy()
        logging("Removing Unprofiled cells from training data", resultdir, "log.txt")
        train_df = train_df.query('genotype_merged != "Unprofiled"')
        logging("Train set Labels before balancing by label:", resultdir, "log.txt")
        logging(train_df["Y"].value_counts(), resultdir, "log.txt")

        # Balance by Label (WT, MUT, False Positive)
        train_df = ml.balance_by_label(
            train_df,
            "Y",
            balance_values=[0, 2],
            min_freq=100,
            random_seed=random_seed,
            visualize=False,
        )
        logging("Train set Labels after balancing by label:", resultdir, "log.txt")
        logging(train_df["Y"].value_counts(), resultdir, "log.txt")
        train_df.to_csv(os.path.join(resultdir, f"{target}_2_train.csv"), index=False)
        test_df.to_csv(os.path.join(resultdir, f"{target}_2_test.csv"), index=False)
        train_df = train_df.drop(columns=columns_to_drop)
        test_df = test_df.drop(columns=columns_to_drop)
        (
            test_res_2,
            fitted_models,
            _best_model_name,
            _y_pred_best,
            _y_pred_prob_best,
            _best_clf,
        ) = mcl.train_classification_model(
            train_df,
            test_df,
            resultdir,
            target_column="Y",
            y_value_of_interest=2,
            models_to_run=models_to_run,
            ensemble_hard=True,
            ensemble_soft=True,
            random_seed=random_seed,
        )
        test_res = pd.concat([test_res_1, test_res_2])
        test_res = test_res[test_res.columns[~test_res.columns.str.contains("best")]]

        # * ================================================ * #
        # * === Gather all predictions and probabilities === * #
        # * ================================================ * #
        pred_cols = test_res.columns[test_res.columns.str.contains("_pred")].tolist()
        prob_cols = test_res.columns[test_res.columns.str.contains("_prob_")].tolist()
        conf_cols = test_res.columns[test_res.columns.str.contains("_conf")].tolist()
        assert test_res.index.isin(target_df.index).all()
        test_res = pd.merge(
            total_df.loc[test_res.index, columns_to_drop],
            test_res[pred_cols + prob_cols + conf_cols + ["Y"]],
            left_index=True,
            right_index=True,
        ).reset_index(drop=True)

        # Label mapping
        if int2label is not None:
            for col in pred_cols:
                test_res[col] = test_res[col].map(int2label)
            test_res["Y"] = test_res["Y"].map(int2label)
            int2label_proba = {
                c: "_".join(c.split("_")[:-1] + [int2label[int(c.split("_")[-1])]])
                for c in prob_cols
            }
            test_res = test_res.rename(columns=int2label_proba)
        test_res.to_csv(
            os.path.join(resultdir, f"{target}_test_results.csv"), index=False
        )

        # * ====================================================== * #
        # * === Adjust predictions based on confidence and FPR === * #
        # * ====================================================== * #
        # ! Requires `expected` column (containing: 0 or 1)
        if test_res["expected"].nunique() == 1:
            print(
                "All cells are either in the expected cell group or all cells are in the non-expected cell group.\nThe prediction cannot be adjusted."
            )
            continue

        # For each of the models, adjust the prediction based on the confidence
        # for the given threshold of MUT prediction probability, reassign the prediction between WT and MUT
        models = [c.replace("_pred", "") for c in pred_cols]
        for model_name in models:
            pred_col = f"{model_name}_pred"
            prob_col = (
                f"{model_name}_confidence"
                if model_name == "ensemble_hard"
                else f"{model_name}_prob_MUT"
            )
            false_positive_rate = (
                test_res.query("expected == 0")[pred_col]
                .value_counts(normalize=True)
                .get("MUT")
            )
            true_positive_rate = (
                test_res.query("expected == 1")[pred_col]
                .value_counts(normalize=True)
                .get("MUT")
            )
            if false_positive_rate is None:
                test_res[f"{pred_col}_adj"] = test_res[pred_col]
                continue
            if true_positive_rate is None:
                true_positive_rate = 0
            specificity = true_positive_rate - false_positive_rate

            # If there are many false positives, adjust the prediction reduce it (less MUT)
            if false_positive_rate > args.max_allowed_fpr:
                # Find the threshold that minimizes the false positive rate without losing too much sensitivity
                for mut_threshold in np.linspace(0.1, 1.0, 50):
                    # values that were originally predicted as MUT but have low confidence are reassigned as FalsePositive
                    mut_mask = test_res[prob_col] >= mut_threshold
                    test_res[f"{pred_col}_tmp"] = np.where(
                        (test_res[pred_col] == "MUT") & ~mut_mask,
                        "FalsePositive",
                        test_res[pred_col],
                    )
                    if "MUT" not in test_res[f"{pred_col}_tmp"].unique():
                        break
                    assert test_res[f"{pred_col}_tmp"].value_counts().get(
                        "MUT"
                    ) <= test_res[pred_col].value_counts().get("MUT")

                    # re-examine the false positive rate
                    false_positive_rate_tmp = (
                        test_res.query("expected == 0")[f"{pred_col}_tmp"]
                        .value_counts(normalize=True)
                        .get("MUT")
                    )
                    true_positive_rate_tmp = (
                        test_res.query("expected == 1")[f"{pred_col}_tmp"]
                        .value_counts(normalize=True)
                        .get("MUT")
                    )
                    if (
                        false_positive_rate_tmp is None
                        or true_positive_rate_tmp is None
                    ):
                        break
                    specificity_tmp = true_positive_rate_tmp - false_positive_rate_tmp
                    if (true_positive_rate - true_positive_rate_tmp < 0.2) and (
                        specificity_tmp > specificity
                    ):
                        test_res[f"{pred_col}_adj"] = test_res[f"{pred_col}_tmp"]
                    if false_positive_rate_tmp < args.max_allowed_fpr:
                        break
            # If false positive rate is low, adjust to get more MUT
            elif false_positive_rate < 0.02:
                # Decrease the threshold to increase the number of MUT but keep the false positive rate low
                for mut_threshold in np.linspace(1.0, 0.01, 50):
                    # values that were originally predicted as WT or FalsePositive but have relatively high confidence are reassigned as MUT
                    mut_mask = test_res[prob_col] >= mut_threshold
                    test_res[f"{pred_col}_tmp"] = np.where(
                        (test_res[pred_col] != "MUT") & mut_mask,
                        "MUT",
                        test_res[pred_col],
                    )
                    false_positive_rate_tmp = (
                        test_res.query("expected == 0")[f"{pred_col}_tmp"]
                        .value_counts(normalize=True)
                        .get("MUT")
                    )
                    true_positive_rate_tmp = (
                        test_res.query("expected == 1")[f"{pred_col}_tmp"]
                        .value_counts(normalize=True)
                        .get("MUT")
                    )
                    if true_positive_rate_tmp is None:
                        true_positive_rate_tmp = 0
                    specificity_tmp = true_positive_rate_tmp - false_positive_rate_tmp
                    if (
                        false_positive_rate_tmp < args.max_allowed_fpr
                        and specificity_tmp >= specificity
                    ):
                        test_res[f"{pred_col}_adj"] = test_res[f"{pred_col}_tmp"]
                        assert true_positive_rate_tmp >= true_positive_rate
            if f"{pred_col}_adj" not in test_res.columns:
                test_res[f"{pred_col}_adj"] = test_res[pred_col]
        test_res.drop(
            columns=test_res.columns[test_res.columns.str.contains("_tmp")],
            inplace=True,
        )

        # === Save Best Model Prediction in 'best_pred_adj' column === #
        performance_dict = {}
        for model_name in models:
            performance = custom_performance_metric(
                test_res["Y"],
                test_res[f"{model_name}_pred_adj"],
                y_value_of_interest="MUT",
                alpha=args.alpha,
            )
            performance_dict[model_name] = performance
        models_sorted = [
            k
            for k, v in sorted(
                performance_dict.items(), key=lambda item: item[1], reverse=True
            )
        ]
        best_model = models_sorted[0]
        test_res["best_pred_adj"] = test_res[f"{best_model}_pred_adj"]
        test_res.to_csv(
            os.path.join(resultdir, f"{target}_test_results_adjusted.csv"), index=False
        )

    # === Aggregate Results === #
    agg_pred_df = []
    for target in targets:
        resultdir = os.path.join(args.outdir, target)
        if not os.path.exists(resultdir):
            if target in skipped_targets:
                print(
                    f"{target} was not denoised using ML. Adding IronThrone results as the final genotypes."
                )
                test_res = total_df[total_df["target"] == target]
                test_res["best_pred_adj"] = test_res["genotype_merged"].copy()
            else:
                raise FileNotFoundError(f"{target} does not have a result directory.")
        else:
            test_res = pd.read_csv(
                os.path.join(resultdir, f"{target}_test_results_adjusted.csv")
            )
        agg_pred_df.append(test_res)
    agg_pred_df = pd.concat(agg_pred_df, join="outer").drop_duplicates()
    createFolder(os.path.join(args.outdir, "aggregated"))
    # agg_pred_df.to_csv(
    #     os.path.join(args.outdir, "aggregated", "prediction_results_aggregated.csv"),
    #     index=False,
    # )

    # === Save Best Model Prediction in 'best_pred_adj' column and recover Unprofiled cells === #
    pred_cols = agg_pred_df.columns[agg_pred_df.columns.str.contains("_pred")].tolist()
    pred_cols = [c for c in pred_cols if "_adj" in c and "best" not in c]
    targets = agg_pred_df["target"].unique()
    final_pred_df = []
    for target in targets:
        test_res = agg_pred_df[agg_pred_df["target"] == target]
        performance_dict = {}
        for pred_col in pred_cols:
            test_res[pred_col] = test_res.apply(
                lambda row: "Unprofiled"
                if row["genotype_merged"] == "Unprofiled" and row[pred_col] != "MUT"
                else row[pred_col],
                axis=1,
            )
            test_res[pred_col] = test_res[pred_col].str.replace("FalsePositive", "WT")
            assert test_res["BC"].isin(total_df["BC"]).all()
            test_res["gt"] = test_res["Y"].str.replace("FalsePositive", "WT")
            test_tmp = test_res.query(f'{pred_col} != "Unprofiled"')
            performance = custom_performance_metric(
                test_tmp["gt"],
                test_tmp[pred_col],
                y_value_of_interest="MUT",
                alpha=args.alpha,
            )
            performance_dict[pred_col] = performance
        models_sorted = [
            k
            for k, v in sorted(
                performance_dict.items(), key=lambda item: item[1], reverse=True
            )
        ]
        best_model = models_sorted[0]
        test_res["best_pred_adj"] = test_res[best_model]
        final_pred_df.append(test_res)
        logging(
            f"Target: {target}, Best Model: {best_model}",
            os.path.join(args.outdir, target),
            "log.txt",
        )
    final_pred_df = pd.concat(final_pred_df)
    final_pred_df["genotype_merged"] = final_pred_df["genotype_merged"].str.replace(
        "Ambiguous", "Unprofiled"
    )

    # Post-process
    final_pred_df["best_pred_adj"] = final_pred_df.apply(
        lambda row: "Unprofiled"
        if row["genotype_merged"] == "Unprofiled" and row["best_pred_adj"] != "MUT"
        else row["best_pred_adj"],
        axis=1,
    )
    final_pred_df["best_pred_adj"] = final_pred_df["best_pred_adj"].str.replace(
        "FalsePositive", "WT"
    )
    assert final_pred_df["BC"].isin(total_df["BC"]).all()
    final_pred_df.loc[
        (final_pred_df["expected"] == 1)
        & (final_pred_df["Y"] == "MUT")
        & (final_pred_df["best_pred_adj"] != "MUT"),
        "Y",
    ] = "WT"
    final_pred_df = final_pred_df[
        (final_pred_df["genotype_merged"].isin(["WT", "MUT"]))
        | (final_pred_df["best_pred_adj"] == "MUT")
    ]
    createFolder(os.path.join(args.outdir, "aggregated"))
    final_pred_df.to_csv(
        os.path.join(
            args.outdir,
            "genotype_refined.csv",
        ),
        index=False,
    )
    print(
        f"Final prediction results written to {os.path.join(args.outdir, 'genotype_refined.csv')}"
    )

    # * =================================== * #
    # * ====== Visualize the results ====== * #
    # * =================================== * #
    import matplotlib.pyplot as plt
    import math

    x_var = "expected"
    final_pred_df["expected"] = final_pred_df["expected"].replace(
        {0: "Not-Expected", 1: "Expected"}
    )

    # === Visualize IronThrone Genoyptes === #
    composition_var = "best_pred_adj"

    targets = final_pred_df["target"].unique()

    n_rows = math.ceil(len(targets) / 2)
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 3.8, n_rows * 3.5))
    gs = GridSpec(
        n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.35
    )  # Adjust hspace and wspace here
    axes = [fig.add_subplot(gs[i, j]) for i in range(n_rows) for j in range(n_cols)]
    for i, target in enumerate(targets):
        target_df = final_pred_df.query(f"target == '{target}'").query(
            f'{composition_var} != "Unprofiled"'
        )  # ! Exclude Unprofiled cells
        expected_sample = target_df["expected_sample"].unique()[0]
        assert target_df["target"].nunique() == 1

        _ = mvz.plot_composition(
            target_df,
            x_var,
            composition_var,
            composition_var_order=["WT", "MUT"],  # "Unprofiled",
            title=f"{target} ({expected_sample})",
            # figsize=(4, 6),
            ax=axes[i],
            # outdir=resultdir,
            palette={"MUT": "#001bc0", "WT": "#8ce", "Unprofiled": "lightgrey"},
            annotate_counts=True,
        )
    plt.savefig(
        os.path.join(
            args.outdir, "refined_genotype_proportions_in_expected_cell_group.png"
        ),
        dpi=150,
        bbox_inches="tight",
    )

    # === Visualize IronThrone Genoyptes === #
    composition_var = "genotype_merged"

    targets = final_pred_df["target"].unique()

    n_rows = math.ceil(len(targets) / 2)
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 3.8, n_rows * 3.5))
    gs = GridSpec(
        n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.35
    )  # Adjust hspace and wspace here
    axes = [fig.add_subplot(gs[i, j]) for i in range(n_rows) for j in range(n_cols)]
    for i, target in enumerate(targets):
        target_df = final_pred_df.query(f"target == '{target}'").query(
            f'{composition_var} != "Unprofiled"'
        )  # ! Exclude Unprofiled cells
        expected_sample = target_df["expected_sample"].unique()[0]
        assert target_df["target"].nunique() == 1

        _ = mvz.plot_composition(
            target_df,
            x_var,
            composition_var,
            composition_var_order=["WT", "MUT"],  # "Unprofiled",
            title=f"{target} ({expected_sample})",
            # figsize=(4, 6),
            ax=axes[i],
            # outdir=resultdir,
            palette={"MUT": "#001bc0", "WT": "#8ce", "Unprofiled": "lightgrey"},
            annotate_counts=True,
        )
    plt.savefig(
        os.path.join(
            args.outdir, "initial_genotype_proportions_in_expected_cell_group.png"
        ),
        dpi=150,
        bbox_inches="tight",
    )

    # * === Visualize Adjusted Prediction Results === #
    targets = final_pred_df["target"].unique()

    n_rows = math.ceil(len(targets) / 2)
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 3.8, n_rows * 3.5))
    gs = GridSpec(
        n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.35
    )  # Adjust hspace and wspace here
    axes = [fig.add_subplot(gs[i, j]) for i in range(n_rows) for j in range(n_cols)]
    for i, target in enumerate(targets):
        target_df = final_pred_df.query(f"target == '{target}'")
        expected_sample = target_df["expected_sample"].unique()[0]
        assert target_df["target"].nunique() == 1

        # Input
        x_var = "sample"  # is_non_mut_group  # "sample"
        composition_var = "best_pred_adj"

        # Calculate performance metrics
        tmp = target_df.copy()
        # ! Redefine label ('gt' column)
        tmp["Y"] = tmp["Y"].str.replace("FalsePositive", "WT")
        tmp["gt"] = tmp["Y"].copy()
        tmp.loc[(tmp["expected"] == 1) & (tmp[composition_var] == "MUT"), "gt"] = "MUT"

        tmp[composition_var] = tmp[composition_var].replace("Unprofiled", "WT")
        acc, prec, recall, f1, conf_matrix = evaluate_prediction(
            tmp,
            pred_col=composition_var,
            y_col="gt",
            labels=["WT", "MUT"],
            pos_label="MUT",
            report=False,
        )
        metrics_report = (
            f"Acc: {acc:.2f},  Prec: {prec:.2f},  Rec: {recall:.2f}, F1: {f1:.2f}"
        )

        target_df = target_df.query(
            f'{composition_var} != "Unprofiled"'
        )  # ! Exclude Unprofiled cells

        _ = mc.plot_composition(
            target_df,
            x_var,
            composition_var,
            composition_var_order=["WT", "MUT"],  # "Unprofiled",
            title=f"{target} ({expected_sample})",
            ax=axes[i],
            palette={"MUT": "#001bc0", "WT": "#8ce", "Unprofiled": "lightgrey"},
            annotate_counts=True,
        )
        x_mid = (axes[i].get_xlim()[0] + axes[i].get_xlim()[1]) / 2
        axes[i].text(
            x=x_mid,
            y=1.075,
            s=metrics_report,
            ha="center",
            size=9,
            # weight="bold",
            color="0.2",
        )
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.savefig(
        os.path.join(
            args.outdir,
            "refined_genotype_proportions_in_expected_cell_group_2.pdf",
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


if __name__ == "__main__":
    main()
