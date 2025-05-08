"""
Modified version of mwpak functions for specific use cases.
Version: v0.1
Modified on: Jan 25, 2025
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
from matplotlib import rcParams
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    f1_score,
    precision_score,
    accuracy_score,
    recall_score,
)


rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.labelsize"] = "7"
rcParams["axes.titlesize"] = "7"
rcParams["font.size"] = "7"
rcParams["xtick.labelsize"] = "7"
rcParams["ytick.labelsize"] = "7"
rcParams["axes.linewidth"] = "0.5"
rcParams["xtick.major.width"] = "0.5"
rcParams["ytick.major.width"] = "0.5"
rcParams["xtick.major.size"] = "3"
rcParams["ytick.major.size"] = "3"
rcParams["mathtext.fontset"] = "custom"
rcParams["mathtext.rm"] = "serif"
rcParams["pdf.fonttype"] = 42
rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"


# ====== Composition Plot (Title position adjusted upward) ====== #
def plot_composition(
    df,
    x_var,
    composition_var,
    composition_var_order=None,
    title=None,
    figsize=(3.5, 3.5),  # Changed default figsize to make room for legend
    outdir=None,
    ax=None,
    top_right_frame=False,
    fontsize=11,
    linewidth=2.0,
    annotate_counts=True,
    **kwargs,
):
    """
    Plot the composition of a categorical variable stratified by another categorical variable.
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the data.
    x_var : str
        The column name of the x variable.
    composition_var : str
        The column name of the composition variable.
    composition_var_order : list, optional
        The order of the composition variable values. The default is None.
    title : str, optional
        The title of the plot. The default is None.
    figsize : tuple, optional
        The figure size. The default is (8, 6).
    outdir : str, optional
        The output directory to save the plot. The default is None.
    ax : matplotlib.axes.Axes, optional
        The axes to plot on. The default is None.
    top_right_frame : bool, optional
        Whether to show the top and right frame. The default is False.
    fontsize : int, optional
        The font size of the labels. The default is 11.
    Returns
    -------
    None.
    """
    # Make sure x_var is a categorical variable
    df[x_var] = df[x_var].astype(str)

    # Get the counts of each x_var value and composition_var value
    df_counts = (
        df.groupby([x_var, composition_var], observed=False)
        .size()
        .unstack(fill_value=0)
    )
    # Calculate the proportions of the composition_var values for each x_var value
    df_proportions = df_counts.div(df_counts.sum(axis=1), axis=0)
    # Convert to long format for seaborn
    df_proportions_long = df_proportions.reset_index().melt(
        id_vars=x_var, var_name=composition_var, value_name="proportion"
    )

    if composition_var_order is not None:
        # Sort df_proportions_long by x_var and composition_var_order
        df_proportions_long[composition_var] = pd.Categorical(
            df_proportions_long[composition_var],
            categories=composition_var_order,
            ordered=True,
        )
        df_proportions_long = df_proportions_long.sort_values(
            by=[x_var, composition_var]
        )
    else:
        # Sort df_proportions_long by x_var to ensure consistent order
        df_proportions_long = df_proportions_long.sort_values(
            by=[x_var, composition_var], ascending=[True, True]
        )

    df_proportions = df_proportions.loc[
        df_proportions_long[x_var].unique(),
        df_proportions_long[composition_var].unique()[::-1],
    ]
    df_counts = df_counts.loc[
        df_proportions_long[x_var].unique(),
        df_proportions_long[composition_var].unique()[::-1],
    ]

    # Create the stacked bar plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Set default palette if not provided in kwargs
    if "palette" not in kwargs:
        kwargs["palette"] = "tab10"

    sns.histplot(
        data=df_proportions_long,
        x=x_var,
        weights="proportion",
        hue=composition_var,
        multiple="stack",
        shrink=0.8,
        lw=linewidth,
        edgecolor="0.2",
        ax=ax,
        **kwargs,
    )

    # Legend formatting
    legend = ax.get_legend()
    if legend is not None:
        labels = [t.get_text() for t in legend.get_texts()]
        if isinstance(kwargs["palette"], str):
            palette = plt.get_cmap(kwargs["palette"])
            colors = [to_hex(palette(i)) for i in range(palette.N)]
        else:
            colors = [kwargs["palette"][l] for l in labels]

        custom_handles = []
        for label, color in zip(labels, colors):
            patch = mpatches.Patch(color=color, label=label)
            custom_handles.append(patch)

        ax.legend(
            handles=custom_handles,
            loc=10,
            bbox_to_anchor=(0.5, -0.15) if len(labels) <= 3 else (0.5, -0.21),
            ncol=3,
            frameon=False,
            fontsize=14,
            labelcolor="0.2",
            prop={"weight": "bold"},
        )

    plt.tight_layout()

    # Bar annotations (with raw values)
    if annotate_counts:
        # Iterate over the DataFrame and plot text annotations
        for i in range(len(df_proportions)):
            # Get the index (x_var, composition_var) for each bar
            idx = df_proportions.index[i]
            # Get the actual raw values from the original data (df_counts)
            raw_values = df_counts.loc[idx].values
            # Bottom position for stacking
            bottom = 0
            # Annotate for each Pred value (0 and 1)
            for j, value in enumerate(raw_values):
                # Calculate position to place the annotation
                text = ax.text(
                    i,
                    bottom + df_proportions.iloc[i, j] / 2,  # Halfway up each bar
                    f"{value:.0f}",
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    weight="bold",
                    color="0.2",
                )
                text.set_path_effects(
                    [
                        PathEffects.Stroke(linewidth=3, foreground="white"),
                        PathEffects.Normal(),
                    ]
                )
                bottom += df_proportions.iloc[
                    i, j
                ]  # Update bottom position for next stack

    # Plot margins
    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(linewidth)
        ax.spines[axis].set_color("0.2")

    if not top_right_frame:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Axis ticks
    ax.tick_params(width=linewidth, color="0.2")
    # Set font size and weight for y-tick labels
    for label in ax.get_yticklabels():
        label.set_fontsize(fontsize)
        label.set_weight("bold")
        label.set_color("0.2")

    # Set font size and weight for x-tick labels
    for label in ax.get_xticklabels():
        label.set_fontsize(fontsize)
        label.set_weight("bold")
        label.set_color("0.2")

    # Title and axis labels
    if title is not None:
        ax.set_title(title, size=fontsize + 2, weight="bold", color="0.2", pad=15)
    else:
        ax.set_title(
            f"{composition_var} Compositions for Each {x_var}",
            size=fontsize + 2,
            weight="bold",
            color="0.2",
        )
    ax.set_xlabel("")  # f"{x_var}"
    ax.set_ylabel(
        f"{composition_var} (%)",
        size=fontsize,
        weight="bold",
        color="0.2",  # - Percentage (%)
        labelpad=-1,
    )
    if outdir is not None:
        fig.savefig(
            os.path.join(
                outdir, f"{composition_var}_compositions_for_each_{x_var}.png"
            ),
            dpi=300,
            bbox_inches="tight",
        )
    return df_counts, df_proportions


def evaluate_prediction(
    pred_df, pred_col="pred", y_col="Y", labels=[0, 1], pos_label=1, report=True
):
    # Validate the labels
    if (
        not pd.Series(pred_df[y_col].unique()).isin(labels).all()
        or not pd.Series(pred_df[pred_col].unique()).isin(labels).all()
        # not pd.Series(labels).isin(pred_df[y_col].unique()).all()
        # or not pd.Series(labels).isin(pred_df[pred_col].unique()).all()
    ):
        print("Labels mismatch between the prediction and the target.")
        print(f"\tLabels Specified: {labels}")
        print(f"\tPrediction labels: {pred_df[pred_col].unique()}")
        print(f"\tTarget labels: {pred_df[y_col].unique()}")
        # acc = prec = recall = f1 = 0
        # return acc, prec, recall, f1, None

    # === Report === #
    # Confusion matrix
    conf_matrix = confusion_matrix(pred_df[y_col], pred_df[pred_col], labels=labels)

    # Classification report
    if report:
        print(classification_report(pred_df[y_col], pred_df[pred_col]))

    # === Metrics scores === #
    # TN, FP, FN, TP
    tn, fp, fn, tp = confusion_matrix(
        pred_df[y_col], pred_df[pred_col], labels=labels
    ).ravel()
    if report:
        print(f"TN: {tn}, FP: {fp}, FN: {fn}, TP: {tp}")

    # Accuracy
    acc_man = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) != 0 else np.nan
    acc = accuracy_score(pred_df[y_col], pred_df[pred_col])
    acc = acc if not np.isnan(acc) else acc_man
    acc_man = acc_man if not np.isnan(acc_man) else acc

    # Precision
    prec_man = tp / (tp + fp) if (tp + fp) != 0 else np.nan
    prec = precision_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    prec = prec if not np.isnan(prec) else prec_man
    prec_man = prec_man if not np.isnan(prec_man) else prec

    # Recall
    recall_man = tp / (tp + fn) if (tp + fn) != 0 else np.nan
    recall = recall_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    recall = recall if not np.isnan(recall) else recall_man
    recall_man = recall_man if not np.isnan(recall_man) else recall

    # F1 score
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
