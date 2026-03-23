import os
import sys
import numpy as np
import pandas as pd
import warnings
import argparse
import matplotlib.pyplot as plt
import math
from matplotlib.gridspec import GridSpec
from sklearn.model_selection import train_test_split

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import utils as ml
import utils.classification as mcl
import utils.vz as mvz
from utils.classification import custom_performance_metric
from utils import createFolder, logging


def format_expected_sample_label(series: pd.Series) -> str:
    expected_samples = sorted(
        {
            value.strip()
            for entry in series.dropna().astype(str)
            for value in entry.split(",")
            if value.strip()
        }
    )
    return ", ".join(expected_samples)

parser = argparse.ArgumentParser(description="ML denoise and genotype adjustment")
parser.add_argument(
    "--input-features",
    dest="input_features",
    type=str,
    default="ironthrone_out_pp_ml_input.csv",
    help="Path to the input features CSV",
)
parser.add_argument(
    "--outdir",
    type=str,
    default="Project/ML_genotyping",
    help="Directory to write outputs",
)
parser.add_argument(
    "--max-allowed-fpr",
    dest="max_allowed_fpr",
    type=float,
    default=0.03,
    help="Maximum allowed false positive rate when adjusting predictions",
)
parser.add_argument(
    "--alpha",
    type=float,
    default=0.8,
    help="Alpha parameter for custom performance metric (higher favors precision)",
)
parser.add_argument(
    "--random-seed",
    dest="random_seed",
    type=int,
    default=1024,
    help="Random seed for reproducibility",
)
args, _ = parser.parse_known_args()

random_seed = args.random_seed  # 1204
np.random.seed(random_seed)
warnings.filterwarnings("ignore")

createFolder(args.outdir)

# ! Features to drop when training the model
columns_to_drop = [
    "BC",
    "target",
    "sample",
    "expected_sample",
    "expected",
    "experiment",
    "is_non_mut",
    "genotype_geneseq",
    "genotype_probebc",
    "genotype_merged",
    "cell_group",
]
# ! Models to run
models_to_run = [
    "logistic_regression",
    "random_forest",
    "knn",
    "naive_bayes",
    "xgboost",
    "mlp",
    "gradient_boosting",
    "adaboost",
    # "hist_gradient_boosting",
    # "svm",
    # "qda",
]
# ! To convert integer labels to back to original labels (if available, default is None)
int2label = {0: "FalsePositive", 1: "WT", 2: "MUT"}

# %% [load-data]
total_df = pd.read_csv(args.input_features)
assert (
    total_df.query('genotype_merged == "MUT"')["Y_num"]
    .value_counts()
    .index.isin([0, 2])
    .all()
)
total_df.drop("Y", axis=1, inplace=True)
total_df.rename(columns={"Y_num": "Y"}, inplace=True)

# %% [targets]
# === Predict for each target (across all experiments) === #
targets = total_df["target"].unique()
print(f"#Targets: {len(targets)}")
print(targets)

# %% [per-target-training]
skipped_targets = []
for target in targets:
    target_df = total_df[total_df["target"] == target]
    non_target_df = total_df[total_df["target"] != target]

    if (target_df["expected"] == 1).sum() <= 1:
        print(
            f"Skipping {target} due to lack of cells (< 2) in expected sample/cell group"
        )
        skipped_targets.append(target)
        continue
    if (target_df["Y"] == 2).sum() <= 1:
        print(
            f"Skipping {target} due to lack of MUT (< 2) in expected sample/cell group"
        )
        skipped_targets.append(target)
        continue
    # skip if there is no false positive mutant cells
    if (target_df["Y"] == 0).sum() < 1:
        print(f"Skipping {target} due to lack of false positive mutants")
        skipped_targets.append(target)
        continue

    resultdir = os.path.join(args.outdir, target)
    createFolder(resultdir)
    logging(f"# === Target: {target} === #", resultdir, "log.txt")
    logging(
        total_df.query("target == @target")["genotype_merged"].value_counts(),
        resultdir,
        "log.txt",
    )

    # * ================================================================================ * #
    # * === First Run (Train: Half of target + Other targets & Test: The other half) === * #
    # * ================================================================================ * #
    # if any of the values other than 2 in 'Y' is equal to 1, duplicate that row to have two rows
    if (target_df["Y"] == 0).sum() == 1:
        target_df = pd.concat([target_df, target_df[target_df["Y"] == 0]])
    if (target_df["Y"] == 1).sum() == 1:
        target_df = pd.concat([target_df, target_df[target_df["Y"] == 1]])

    target_df_1, target_df_2 = train_test_split(
        target_df,
        test_size=0.5,
        stratify=target_df["Y"],
        random_state=random_seed,
    )

    logging(
        f"First Run: Train: {target_df_1.shape}, Test: {target_df_2.shape}",
        resultdir,
        "log.txt",
    )

    # Train: Half of target + Other targets, Test: The other half
    train_df = pd.concat([target_df_1, non_target_df])
    test_df = target_df_2.copy()

    # ! Removing Unprofiled cells from Training data
    logging("Removing Unprofiled cells from training data", resultdir, "log.txt")
    train_df = train_df.query('genotype_merged != "Unprofiled"')
    test_df = test_df.query('genotype_merged != "Unprofiled"')  # !
    logging("Train set Labels before balancing by label:", resultdir, "log.txt")
    logging(train_df["Y"].value_counts(), resultdir, "log.txt")

    # * Balance by Label (WT, MUT, False Positive)
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

    # === Drop columns that are not features === #
    columns_to_drop = [c for c in columns_to_drop if c in train_df.columns]
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

    # ! Removing Unprofiled cells from Training data
    logging("Removing Unprofiled cells from training data", resultdir, "log.txt")
    train_df = train_df.query('genotype_merged != "Unprofiled"')
    test_df = test_df.query('genotype_merged != "Unprofiled"')  #!
    logging("Train set Labels before balancing by label:", resultdir, "log.txt")
    logging(train_df["Y"].value_counts(), resultdir, "log.txt")

    # * Balance by Label (WT, MUT, False Positive)
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

    # === Drop columns that are not features === #
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
    # ! Remove 'best' columns
    test_res = test_res[test_res.columns[~test_res.columns.str.contains("best")]]

    # * ================================================ * #
    # * === Gather all predictions and probabilities === * #
    # * ================================================ * #
    pred_cols = test_res.columns[test_res.columns.str.contains("_pred")].tolist()
    prob_cols = test_res.columns[test_res.columns.str.contains("_prob_")].tolist()
    conf_cols = test_res.columns[test_res.columns.str.contains("_conf")].tolist()

    assert test_res.index.isin(target_df.index).all()
    test_res = pd.merge(
        total_df.loc[
            test_res.index,
            columns_to_drop,
        ],
        test_res[pred_cols + prob_cols + conf_cols + ["Y"]],
        left_index=True,
        right_index=True,
    ).reset_index(drop=True)

    # * === Label mapping === #
    if int2label is not None:
        # Map integer to label for each prediction column
        for col in pred_cols:
            test_res[col] = test_res[col].map(int2label)
        test_res["Y"] = test_res["Y"].map(int2label)

        # Map integer to label in the pred_proba column names
        int2label_proba = {
            c: "_".join(c.split("_")[:-1] + [int2label[int(c.split("_")[-1])]])
            for c in prob_cols
        }
        test_res = test_res.rename(columns=int2label_proba)

    test_res.to_csv(os.path.join(resultdir, f"{target}_test_results.csv"), index=False)

    # * ================================================================= * #
    # * ====== Adjust the prediction based on the confidence score ====== * #
    # * ================================================================= * #
    # For each of the models, adjust the prediction based on the confidence
    # for the given threshold of MUT prediction probability, reassign the prediction between WT and MUT
    models = [c.replace("_pred", "") for c in pred_cols]

    # ! Requires `expected` column (containing: 0 or 1)
    # if all cells are in the expected cell group, or all cells are in the non-expected cell group,
    # then the prediction cannot be adjusted
    if test_res["expected"].nunique() == 1:
        print(
            "All cells are either in the expected cell group or all cells are in the non-expected cell group.\n"
            "The prediction cannot be adjusted."
        )
        # Just save the original predictions as adjusted predictions
        for model_name in models:
            pred_col = f"{model_name}_pred"
            if f"{pred_col}_adj" not in test_res.columns:
                test_res[f"{pred_col}_adj"] = test_res[pred_col]

    else:
        for model_name in models:
            pred_col = f"{model_name}_pred"
            if model_name == "ensemble_hard":
                prob_col = f"{model_name}_confidence"
            else:
                prob_col = f"{model_name}_prob_MUT"

            # Calculate false positive rate
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

            # * If there are many false positives, adjust the prediction reduce it (less MUT)
            if false_positive_rate > args.max_allowed_fpr:
                # Find the threshold that minimizes the false positive rate without losing too much sensitivity
                for mut_threshold in np.linspace(
                    0.1, 1.0, 50
                ):  # stricter (less MUT: Some MUT -> FalsePositive)
                    mut_mask = test_res[prob_col] >= mut_threshold

                    # values that were originally predicted as MUT but have low confidence are reassigned as FalsePositive
                    test_res[f"{pred_col}_tmp"] = np.where(
                        (test_res[pred_col] == "MUT") & ~mut_mask,
                        "FalsePositive",
                        test_res[pred_col],
                    )
                    if "MUT" not in test_res[f"{pred_col}_tmp"].unique():
                        break
                    assert test_res[f"{pred_col}_tmp"].value_counts().get(
                        "MUT"
                    ) <= test_res[f"{pred_col}"].value_counts().get("MUT")

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

                    # Only adjust if we didn't lose too much true positive MUT calls and specificity is improved
                    if (true_positive_rate - true_positive_rate_tmp < 0.2) and (
                        specificity_tmp > specificity
                    ):
                        test_res[f"{pred_col}_adj"] = test_res[f"{pred_col}_tmp"]
                    if false_positive_rate_tmp < args.max_allowed_fpr:
                        break
            # * If false positive rate is low, adjust to get more MUT
            elif false_positive_rate < 0.02:
                # Decrease the threshold to increase the number of MUT but keep the false positive rate low
                for mut_threshold in np.linspace(
                    1.0, 0.01, 50
                ):  # looser (more MUT: Some WT/FaslePositive -> MUT)
                    mut_mask = test_res[prob_col] >= mut_threshold

                    # values that were originally predicted as WT or FalsePositive but have relatively high confidence are reassigned as MUT
                    test_res[f"{pred_col}_tmp"] = np.where(
                        (test_res[pred_col] != "MUT") & mut_mask,
                        "MUT",
                        test_res[pred_col],
                    )

                    # examine false & true positive rates
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

                    # Only adjust false positive is still low and specificity is improved
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

    # * === Save Best Model Prediction in 'genotype_final' column === #
    # Sort models by performance
    performance_dict = {}
    for model_name in models:
        performance = custom_performance_metric(
            test_res["Y"],
            test_res[
                f"{model_name}_pred_adj"
            ],  # ! For sorting, use the "adjusted" prediction
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
    test_res["genotype_final"] = test_res[f"{best_model}_pred_adj"]

    test_res.to_csv(
        os.path.join(resultdir, f"{target}_test_results_adjusted.csv"), index=False
    )

# %% [aggregate]
# * =============================== * #
# * ====== Aggregate Results ====== * #
# * =============================== * #
agg_pred_df = []
for target in targets:
    resultdir = os.path.join(args.outdir, target)

    # check if resultdir exists
    if not os.path.exists(resultdir):
        if target in skipped_targets:
            print(
                f"{target} was not denoised using ML. Adding IronThrone results as the final genotypes."
            )
            test_res = total_df[total_df["target"] == target]
            test_res["genotype_final"] = test_res["genotype_merged"].copy()
            test_res["Y"] = test_res["Y"].map(int2label)
        else:
            raise FileNotFoundError(f"{target} does not have a result directory.")
    else:
        test_res = pd.read_csv(
            os.path.join(resultdir, f"{target}_test_results_adjusted.csv")
        )
    agg_pred_df.append(test_res)

agg_pred_df = pd.concat(agg_pred_df, join="outer").drop_duplicates()

createFolder(os.path.join(args.outdir, "aggregated"))
agg_pred_df.to_csv(
    os.path.join(args.outdir, "aggregated", "prediction_results_aggregated.csv"),
    index=False,
)

# * ====== Find & Save Best Model ====== * #
pred_cols = agg_pred_df.columns[agg_pred_df.columns.str.contains("_pred")].tolist()
pred_cols = [c for c in pred_cols if "_adj" in c and "best" not in c]
targets = agg_pred_df["target"].unique()

final_pred_df = []
for target in targets:
    test_res = agg_pred_df[agg_pred_df["target"] == target]

    # If target was not refined, just use IronThrone results
    if target in skipped_targets:
        final_pred_df.append(test_res)
        print(f"Target: {target}, Best Model: IronThrone Result (Skipped)")
        continue

    performance_dict = {}
    for pred_col in pred_cols:
        test_res[pred_col] = test_res.apply(
            lambda row: "Unprofiled"
            if row["genotype_merged"] == "Unprofiled" and row[pred_col] != "MUT"
            else row[pred_col],
            axis=1,
        )
        # ! Keep Original pred
        test_res[f"{pred_col}_tmp"] = test_res[pred_col].str.replace(
            "FalsePositive", "WT"
        )
        test_res["gt"] = test_res["Y"].str.replace("FalsePositive", "WT")
        assert test_res["BC"].isin(total_df["BC"]).all()
        test_tmp = test_res.query(f'{pred_col}_tmp != "Unprofiled"')

        performance = custom_performance_metric(
            test_tmp["gt"],
            test_tmp[pred_col],
            y_value_of_interest="MUT",
            alpha=args.alpha,
        )
        performance_dict[pred_col] = performance

    # Sort models by performance
    models_sorted = [
        k
        for k, v in sorted(
            performance_dict.items(), key=lambda item: item[1], reverse=True
        )
    ]
    # Save Best Model Prediction in 'genotype_final' column
    best_model = models_sorted[0]
    test_res["genotype_final"] = test_res[best_model]
    final_pred_df.append(test_res)
    print(f"Target: {target}, Best Model: {best_model}")
final_pred_df = pd.concat(final_pred_df)
# final_pred_df["genotype_merged"] = final_pred_df["genotype_merged"].str.replace(
#     "Ambiguous", "Unprofiled"
# )

# === Drop Ambiguous/Unprofiled cells === #
final_pred_df = final_pred_df[final_pred_df["genotype_merged"] != "Ambiguous"]
final_pred_df = final_pred_df[final_pred_df["genotype_final"] != "Unprofiled"]

# === Drop temporary columns === #
# drop columns with '_tmp' in the column names
final_pred_df = final_pred_df.loc[:, ~final_pred_df.columns.str.contains("_tmp")]
# drop 'gt' column
final_pred_df = final_pred_df.drop(columns=["gt"], errors="ignore")


# %% [post-process]
# * ====== Post-Process ====== * #
assert (final_pred_df["genotype_merged"] != "Unprofiled").all()
assert final_pred_df["BC"].isin(total_df["BC"]).all()

# Y values of the cells that were called MUT in the expected cell group by ironthrone but are predicted as WT by the model are set to FalsePositive
# (Since they are false positive noises and saved separately as 'gt')
final_pred_df["gt"] = final_pred_df["Y"].copy()
final_pred_df.loc[
    (final_pred_df["expected"] == 1)
    & (final_pred_df["Y"] == "MUT")
    & (final_pred_df["genotype_final"] != "MUT"),
    "gt",
] = "FalsePositive"

createFolder(os.path.join(args.outdir, "aggregated"))
final_output_path = os.path.join(args.outdir, "got_multi_out_refined.csv")
final_pred_df.to_csv(final_output_path, index=False)
print(f"Final prediction results written to {final_output_path}")

# %% [viz-refined]
# * =================================== * #
# * ====== Visualize the results ====== * #
# * =================================== * #
x_var = "expected"
final_pred_df["expected"] = final_pred_df["expected"].replace(
    {0: "Not-expected", 1: "Expected"}
)

# === Visualize Adjusted Genotypes === #
composition_var = "genotype_final"

targets = final_pred_df["target"].unique()

n_rows = math.ceil(len(targets) / 2) if len(targets) > 1 else 1
n_cols = 2 if len(targets) > 1 else 1
fig = plt.figure(figsize=(n_cols * 3.8, n_rows * 3.5))
gs = GridSpec(
    n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.35
)  # Adjust hspace and wspace here
axes = [fig.add_subplot(gs[i, j]) for i in range(n_rows) for j in range(n_cols)]
for i, target in enumerate(targets):
    target_df = final_pred_df.query(f"target == '{target}'").query(
        f'{composition_var} != "Unprofiled"'
    )  # ! Exclude Unprofiled cells
    expected_sample = format_expected_sample_label(target_df["expected_sample"])
    assert target_df["target"].nunique() == 1

    _ = mvz.plot_composition(
        target_df,
        x_var,
        composition_var,
        composition_var_order=["WT", "FalsePositive", "MUT"],  # "Unprofiled",
        title=f"{target} ({expected_sample})",
        # figsize=(4, 6),
        ax=axes[i],
        # outdir=resultdir,
        palette={"MUT": "#001bc0", "WT": "#8ce", "FalsePositive": "lightgrey"},
        annotate_counts=True,
    )
plt.savefig(
    os.path.join(
        args.outdir, "refined_genotype_proportions_in_expected_cell_group.pdf"
    ),
    dpi=150,
    bbox_inches="tight",
)

# %% [viz-initial]
# === Visualize IronThrone Genotypes === #
composition_var = "genotype_merged"
x_var = "expected"

targets = final_pred_df["target"].unique()

n_rows = math.ceil(len(targets) / 2) if len(targets) > 1 else 1
n_cols = 2 if len(targets) > 1 else 1
fig = plt.figure(figsize=(n_cols * 3.8, n_rows * 3.5))
gs = GridSpec(
    n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.35
)  # Adjust hspace and wspace here
axes = [fig.add_subplot(gs[i, j]) for i in range(n_rows) for j in range(n_cols)]
for i, target in enumerate(targets):
    target_df = final_pred_df.query(f"target == '{target}'").query(
        f'{composition_var} != "Unprofiled"'
    )  # ! Exclude Unprofiled cells
    # target_df = target_df.query(
    #     f'{composition_var} != "FalsePositive"'
    # )  # ! Exclude FalsePositive cells
    expected_sample = format_expected_sample_label(target_df["expected_sample"])
    assert target_df["target"].nunique() == 1

    _ = mvz.plot_composition(
        target_df,
        x_var,
        composition_var,
        composition_var_order=["FalsePositive", "WT", "MUT"],  # "Unprofiled",
        title=f"{target} ({expected_sample})",
        # figsize=(4, 6),
        ax=axes[i],
        # outdir=resultdir,
        palette={"MUT": "#001bc0", "WT": "#8ce", "Unprofiled": "lightgrey"},
        annotate_counts=True,
    )
plt.savefig(
    os.path.join(
        args.outdir, "initial_genotype_proportions_in_expected_cell_group.pdf"
    ),
    dpi=150,
    bbox_inches="tight",
)
