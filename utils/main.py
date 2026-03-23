import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)

warnings.filterwarnings("ignore")


def logging(msg, outdir, log_fpath):
    fpath = os.path.join(outdir, log_fpath)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(fpath, "a") as fw:
        fw.write("%s\n" % msg)
    print(msg)


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print("Error: Creating directory. " + directory)


def balance_by_label(
    train_df,
    balance_column="Y",
    balance_values=None,
    random_seed=42,
    min_freq=100,
    visualize=False,
    outdir=None,
):
    """Balance a dataframe by oversampling or undersampling the requested label column."""
    if balance_values is None:
        max_count = train_df[balance_column].value_counts().max()
    else:
        max_count = (
            train_df[train_df[balance_column].isin(balance_values)][balance_column]
            .value_counts()
            .max()
        )

    if max_count < min_freq:
        max_count = min_freq

    balance_values = train_df[balance_column].unique()

    balanced_dfs = []
    for value in balance_values:
        value_df = train_df[train_df[balance_column] == value]
        if value_df.shape[0] < max_count:
            value_df = value_df.sample(
                n=max_count, replace=True, random_state=random_seed
            )
        elif value_df.shape[0] > max_count:
            value_df = value_df.sample(
                n=max_count, replace=False, random_state=random_seed
            )
        balanced_dfs.append(value_df)

    balanced_df = pd.concat(balanced_dfs)
    assert balanced_df[balance_column].value_counts().nunique() == 1

    if visualize:
        plt.figure(figsize=(8, 6))
        sns.countplot(data=balanced_df, x=balance_column, order=balance_values)
        plt.title(f"Balanced Data Composition by {balance_column}")
        plt.xlabel(balance_column)
        plt.ylabel("Count")
        if outdir:
            plt.savefig(
                os.path.join(outdir, f"balanced_data_composition_{balance_column}.png"),
                dpi=300,
                bbox_inches="tight",
            )
        plt.show()

    return balanced_df


def evaluate_prediction(
    pred_df, pred_col="pred", y_col="Y", labels=[0, 1], pos_label=1, report=True
):
    if (
        not pd.Series(pred_df[y_col].unique()).isin(labels).all()
        or not pd.Series(pred_df[pred_col].unique()).isin(labels).all()
    ):
        print("Labels mismatch between the prediction and the target.")
        print(f"\tLabels Specified: {labels}")
        print(f"\tPrediction labels: {pred_df[pred_col].unique()}")
        print(f"\tTarget labels: {pred_df[y_col].unique()}")

    conf_matrix = confusion_matrix(pred_df[y_col], pred_df[pred_col], labels=labels)

    if report:
        print(classification_report(pred_df[y_col], pred_df[pred_col]))

    tn, fp, fn, tp = confusion_matrix(
        pred_df[y_col], pred_df[pred_col], labels=labels
    ).ravel()
    if report:
        print(f"TN: {tn}, FP: {fp}, FN: {fn}, TP: {tp}")

    acc_man = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) != 0 else np.nan
    acc = accuracy_score(pred_df[y_col], pred_df[pred_col])
    acc = acc if not np.isnan(acc) else acc_man
    acc_man = acc_man if not np.isnan(acc_man) else acc

    prec_man = tp / (tp + fp) if (tp + fp) != 0 else np.nan
    prec = precision_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    prec = prec if not np.isnan(prec) else prec_man
    prec_man = prec_man if not np.isnan(prec_man) else prec

    recall_man = tp / (tp + fn) if (tp + fn) != 0 else np.nan
    recall = recall_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    recall = recall if not np.isnan(recall) else recall_man
    recall_man = recall_man if not np.isnan(recall_man) else recall

    f1_man = (
        2 * (prec_man * recall_man) / (prec_man + recall_man)
        if (prec_man + recall_man) != 0
        and not (np.isnan(prec_man) or np.isnan(recall_man))
        else np.nan
    )
    f1 = f1_score(pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label)
    f1 = f1 if not np.isnan(f1) else f1_man
    f1_man = f1_man if not np.isnan(f1_man) else f1

    if report:
        print(f"F1: {f1}, F1_man: {f1_man}")
        print(f"Precision: {prec}, Recall: {recall}")

    return acc, prec, recall, f1, conf_matrix
