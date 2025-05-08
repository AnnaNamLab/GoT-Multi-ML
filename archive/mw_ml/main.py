"""v1.3 - 250402"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


def logging(msg, outdir, log_fpath):
    fpath = os.path.join(outdir, log_fpath)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(fpath, "a") as fw:
        fw.write("%s\n" % msg)
    print(msg)


def balance_by_label(
    train_df,
    balance_column="Y",
    balance_values=None,
    random_seed=42,
    min_freq=100,
    visualize=False,
    outdir=None,
):
    """
    Balance the data by oversampling or undersampling based on the specified column.

    Parameters:
    - train_df (pd.DataFrame): The training dataframe.
    - balance_column (str): The column to balance the data on. Default is 'Y'.
    - balance_values (list): The specific values to consider for calculating the reference number. Default is None.
    - random_seed (int): The random seed for reproducibility. Default is 42.
    - min_freq (int): The minimum frequency threshold for upsampling. Default is 100.
    - visualize (bool): Whether to visualize the resulting data label composition. Default is False.
    - outdir (str): The directory to save the visualization. Default is None.

    Returns:
    - pd.DataFrame: The balanced dataframe.
    """
    # Determine the maximum count for balancing
    if balance_values is None:
        max_count = train_df[balance_column].value_counts().max()
    else:
        max_count = (
            train_df[train_df[balance_column].isin(balance_values)][balance_column]
            .value_counts()
            .max()
        )

    # Ensure the max_count is at least min_freq
    if max_count < min_freq:
        max_count = min_freq

    balance_values = train_df[balance_column].unique()

    balanced_dfs = []
    for value in balance_values:
        value_df = train_df[train_df[balance_column] == value]
        # Upsample or downsample to match max_count
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

    # Visualize the resulting data label composition if requested
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
