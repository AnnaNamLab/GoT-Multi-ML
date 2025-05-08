"""v2.1 - 250309"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
from matplotlib import rcParams

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


# ====== Plotting Composition of a Categorical Variable stratified by another Categorical Variable ====== #
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
    import matplotlib.patheffects as PathEffects
    import matplotlib.patches as mpatches
    from matplotlib.colors import to_hex

    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"

    # Make sure x_var is a categorical variable
    df[x_var] = df[x_var].astype(str)
    df[composition_var] = df[composition_var].astype(str)

    # Get the counts of each x_var value and composition_var value
    df_counts = (
        df.groupby([x_var, composition_var], observed=False)
        .size()
        .unstack(fill_value=0)
    )
    # Calculate the proportions of the composition_var values for each x_var value
    df_proportions = df_counts.div(df_counts.sum(axis=1), axis=0) * 100
    # Convert to long format for seaborn
    df_proportions_long = df_proportions.reset_index().melt(
        id_vars=x_var, var_name=composition_var, value_name="proportion"
    )

    if composition_var_order is not None:
        composition_var_order = [str(x) for x in composition_var_order]
        composition_var_order = [
            x for x in composition_var_order if x in df[composition_var].unique()
        ]

        # Include only the values in composition_var_order
        composition_var_category_order = [
            v
            for v in df_proportions_long[composition_var].unique()
            if v not in composition_var_order
        ] + composition_var_order
        if len(composition_var_category_order) != len(composition_var_order):
            print(
                "Warning: composition_var_order does not contain all the unique values in the data."
            )

        # Sort df_proportions_long by x_var and composition_var_order
        df_proportions_long[composition_var] = pd.Categorical(
            df_proportions_long[composition_var],
            categories=composition_var_category_order,
            ordered=True,
        )
        df_proportions_long = df_proportions_long.sort_values(
            by=[x_var, composition_var]
        )
        df_proportions_long = df_proportions_long[
            df_proportions_long[composition_var].isin(composition_var_order)
        ]
        df_counts = df_counts[composition_var_order]

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
        ax.set_title(title, size=fontsize + 2, weight="bold", color="0.2")
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
        color="0.2",
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


# ====== Plotting Ratio in Stacked Bar Plot with Raw Values Annotations (Ratios, Counts Precalculated) ====== #
def plot_stacked_barplot_precalculated(
    df,
    x="sample",
    y_cols=["mut_ratio", "wt_ratio"],
    annot_cols=["mut_counts", "wt_counts"],
    figsize=(6, 5),
    fontsize=10,
    palette=None,
    legend_labels=None,
    legend_title="",
    title="",
    ax=None,
):
    """
    Create a stacked bar plot showing the distribution of data and the raw counts across a categorical variable.
    Designed for df with ratio and counts already calculated.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the data to visualize with columns for cell states,
        cell type ratios, and cell counts.
    x : str, default="cell_state"
        Column name in df for x-axis categories (cell states).
    y_cols : list, default=["mut_ratio", "wt_ratio"]
        Column names in df for the stacked bar values (proportions).
    annot_cols : list, default=["n_cellstate_mut_cells", "n_cellstate_wt_cells"]
        Column names in df for annotation values (cell counts).
    legend_labels : list, default=None
        Labels to display in the legend. If None, uses the labels from the plot.
    fontsize : int, default=10
        Base font size for text elements.
    palette : str, list, or None, default=None
        Colors for the bars. Can be:
        - None: Uses default seaborn color palette
        - str: Uses seaborn color palette with this name
        - list/iterable: Uses these colors directly
    figsize : tuple, default=(6, 5)
        Figure dimensions (width, height) in inches.
    title : str, default=""
        Title for the plot.
    legend_title : str, default=""
        Title for the legend.

    Returns:
    --------
    ax : axis object
        Axis object for further customization if needed.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as PathEffects
    import seaborn as sns
    import collections

    # Handle palette options
    if palette is None:
        # Use default seaborn palette if None is provided
        palette = sns.color_palette("tab10", len(y_cols))
    elif isinstance(palette, str):
        # If palette is a string, use it as a seaborn palette name
        palette = sns.color_palette(palette, len(y_cols))
    elif not isinstance(palette, str) and isinstance(palette, collections.abc.Iterable):
        # If palette is a non-string iterable, subset it to y_cols length
        palette = list(palette)[: len(y_cols)]

    # Create figure and axis
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    plt.grid(False)

    # Extract the x-axis values
    x_vals = df[x]

    # Plot the stacked bars
    bottom = df[y_cols[0]].copy()
    for i, y_col in enumerate(y_cols):
        y_vals = df[y_col]
        if i == 0:
            ax.bar(
                x_vals,
                y_vals,
                label=y_col,
                color=palette[i],
                edgecolor="0.2",
                lw=2,
                width=0.85,
            )
        else:
            ax.bar(
                x_vals,
                y_vals,
                bottom=bottom,
                label=y_col,
                color=palette[i],
                edgecolor="0.2",
                lw=2,
                width=0.85,
            )
            bottom += y_vals

    # Add annotations for cell counts
    for i, row in df.iterrows():
        cell_state = row[x]

        assert len(y_cols) == len(annot_cols), (
            "y_cols and annot_cols must have the same length"
        )
        bottom = 0
        for j in range(len(y_cols)):
            y_col = y_cols[j]
            annot_col = annot_cols[j]

            text = ax.text(
                x=cell_state,
                y=bottom + (row[y_col] / 2),
                s=f"{row[annot_col]:.0f}",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=fontsize,
                weight="bold",
                color="0.2",
            )
            text.set_path_effects(
                [
                    PathEffects.Stroke(linewidth=2.0, foreground="white"),
                    PathEffects.Normal(),
                ]
            )
            bottom += row[y_col]

    # Set the y-axis limit to 0-1 since we are plotting proportions
    ax.set_ylim(0, 1.001)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(
        handles,
        legend_labels if legend_labels is not None else labels,
        ncol=1,
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        labelcolor="0.2",
        prop={"weight": "bold", "size": fontsize},
        title=legend_title,
        title_fontsize=fontsize,
    )
    leg.get_title().set_weight("bold")

    # Axis ticks
    ax.tick_params(width=1.5, color="0.2")
    plt.xticks(
        size=fontsize,
        rotation=35,
        rotation_mode="anchor",
        ha="right",
        weight="bold",
        color="0.2",
    )
    plt.yticks(size=fontsize, weight="bold", color="0.2")

    # Add labels and title
    plt.xlabel("Cell State", size=fontsize + 2, weight="bold", color="0.2")
    plt.ylabel("Cell Proportion", size=fontsize + 2, weight="bold", color="0.2")

    # Axis styling
    n_bars = len(x_vals)
    plt.xlim(-0.60, n_bars - 0.37)  # Tighten x-axis limits to reduce edge spacing

    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(1.5)
        ax.spines[axis].set_color("0.2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.title(
        title,
        size=fontsize + 3,
        weight="bold",
        color="0.2",
        pad=10,
    )

    plt.tight_layout()

    return ax


def plot_volcano(
    data,
    log2fc="log2FoldChange",
    padj="padj",
    symbol="symbol",
    baseMean=None,
    pval_thresh=0.05,
    log2fc_thresh=0.75,
    to_label=5,
    color_dict=None,
    shape_dict=None,
    fontsize=10,
    colors=["dimgrey", "lightgrey", "black"],
    top_right_frame=False,
    figsize=(5, 5),
    legend_pos=(1.4, 1),
    point_sizes=(15, 150),
    save=False,
    shapes=None,
    shape_order=None,
    **kwargs,
):
    """
    Make a volcano plot from a pandas dataframe of directly from a csv.

    data : pandas.DataFrame or path to csv
    log2fc : string
        column name of log2 Fold-Change values
    padj : string
        column name of the p values to be converted to -log10 P values
    symbol : string
        column name of gene IDs to use
    baseMean : string
        column name of base mean values for each gene. If this is passed,
        the size of the points will vary.
    pval_thresh : numeric
        threshold padj for points to be significant. Also controls horizontal
        line.
    log2fc_thresh : numeric
        threshold for the absolute value of the log2 fold change to be considered
        significant. Also controls vertical lines
    to_label : int or list
        If an int is passed, that number of top down and up genes will be labeled.
        If a list of gene Ids is passed, only those will be labeled.
    color_dict : dictionary
        dictionary to color dots by. Up to 11 categories with default colors.
        Pass list of genes and the category to group them by. {category : ['gene1', gene2]}
        Default colors are: ['dimgrey', 'lightgrey', 'tab:blue', 'tab:orange',
        'tab:green', 'tab:red', 'tab:purple','tab:brown', 'tab:pink',
        'tab:olive', 'tab:cyan']
    shape_dict : dictionary
        dictionary to shape dots by. Up to 6 categories. Pass list of genes as values
        and category as key. {category : ['gene1', gene2], category2 : ['gene3']}
    fontsize : int
        size of labels
    colors : list
        order and colors to use. Default ['dimgrey', 'lightgrey', 'black']
    top_right_frame : Boolean
        Show the top and right frame. True/False
    figsize : tuple
        Size of figure. (x, y)
    point_sizes : tuple
        lower and upper bounds of point sizes. If baseMean is not None.
        (lower, upper)
    save : boolean | string
        If true saves default file name. Pass string as path to output file. Will
        add a .svg/.png to string. Saves as both png and svg.
    shapes : list
        pass matplotlib marker ids to change default shapes/order
        Default shapes order is: ['o', '^', 's', 'X', '*', 'd']
    shape_order : list
        If you want to change the order of categories for your shapes. Pass
        a list of your categories.

    """
    from adjustText import adjust_text

    if isinstance(data, str):
        df = pd.read_csv(data)
    else:
        df = data.copy(deep=True)

    # clean and imput 0s
    df = df.dropna()
    if df[padj].min() == 0:
        print("0s encountered for p value, imputing 1e-323")
        print("impute your own value if you want to avoid this")
        df[padj][df[padj] == 0] = 1e-323

    pval_thresh = -np.log10(pval_thresh)  # convert p value threshold to nlog10
    df["nlog10"] = -np.log10(df[padj])  # make nlog10 column
    df["sorter"] = df["nlog10"] * np.sign(
        df[log2fc]
    )  # compute score: -log10(p_val) * sign(avg_log2FC) - to pick top genes

    # size the dots by basemean if a column id is passed
    if baseMean is not None:
        df["logBaseMean"] = np.log(df[baseMean])
        baseMean = "logBaseMean"
    else:
        point_sizes = None

    # === color dots by {label:[list of genes]} === #
    # make label list of top x genes up and down, or based on list input
    if isinstance(to_label, int):
        label_df = pd.concat(
            (df.sort_values("sorter")[-to_label:], df.sort_values("sorter")[0:to_label])
        )
    else:
        label_df = df[df[symbol].isin(to_label)]

    # color light grey if below thresh, color picked black
    def map_color_simple(a):
        log2FoldChange, zymbol, nlog10 = a
        if zymbol in label_df[symbol].tolist():
            return "picked"

        if abs(log2FoldChange) < log2fc_thresh or nlog10 < pval_thresh:
            return "not DE"
        return "DE"

    if color_dict is None:
        df["color"] = df[[log2fc, symbol, "nlog10"]].apply(map_color_simple, axis=1)
        hues = ["DE", "not DE", "picked"][: len(df.color.unique())]  # order of colors

    # coloring if dictionary passed
    def map_color_complex(a):
        log2FoldChange, zymbol, nlog10 = a

        for k in list(color_dict):
            if zymbol in color_dict[k]:
                return k
        if abs(log2FoldChange) < log2fc_thresh or nlog10 < pval_thresh:
            return "not DE"
        return "DE"

    if color_dict is not None:
        df["color"] = df[[log2fc, symbol, "nlog10"]].apply(map_color_complex, axis=1)
        user_added_cats = sorted(
            [x for x in df.color.unique() if x not in ["DE", "not DE"]]
        )
        hues = ["DE", "not DE"] + user_added_cats
        # hues = hues[: len(df.color.unique())]  # order of colors
        if colors == ["dimgrey", "lightgrey", "black"]:
            colors = [
                "dimgrey",
                "lightgrey",
                "tab:blue",
                "tab:orange",
                "tab:green",
                "tab:red",
                "tab:purple",
                "tab:brown",
                "tab:pink",
                "tab:olive",
                "tab:cyan",
            ]

    # map shapes if dictionary exists
    def map_shape(zymbol):
        for k in list(shape_dict):
            if zymbol in shape_dict[k]:
                return k

        return "other"

    if shape_dict is not None:
        df["shape"] = df[symbol].map(map_shape)
        user_added_cats = [x for x in df["shape"].unique() if x != "other"]
        shape_order = ["other"] + user_added_cats
        if shapes is None:
            shapes = ["o", "^", "s", "X", "*", "d"]
        shapes = shapes[: len(df["shape"].unique())]
        shape_col = "shape"
    else:
        shape_col = None

    # build palette
    colors = colors[: len(hues)]

    # === Plot === #
    plt.figure(figsize=figsize)
    ax = sns.scatterplot(
        data=df,
        x=log2fc,
        y="nlog10",
        hue="color",
        hue_order=hues,
        palette=colors,
        size=baseMean,
        sizes=point_sizes,
        style=shape_col,
        style_order=shape_order,
        markers=shapes,
        **kwargs,
    )

    # make labels
    texts = []
    for i in range(len(label_df)):
        txt = plt.text(
            x=label_df.iloc[i][log2fc],
            y=label_df.iloc[i].nlog10,
            s=label_df.iloc[i][symbol],
            fontsize=fontsize,
            weight="bold",
        )

        txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground="w")])
        texts.append(txt)
    adjust_text(
        texts, expand=(1.1, 1.5), arrowprops=dict(arrowstyle="-", color="k", zorder=5)
    )

    # plot vertical and horizontal lines
    ax.axhline(pval_thresh, zorder=0, c="k", lw=2, ls="--")
    ax.axvline(log2fc_thresh, zorder=0, c="k", lw=2, ls="--")
    ax.axvline(log2fc_thresh * -1, zorder=0, c="k", lw=2, ls="--")

    # make things pretty
    for axis in ["bottom", "left", "top", "right"]:
        ax.spines[axis].set_linewidth(2)

    if not top_right_frame:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax.tick_params(width=2)
    plt.xticks(size=11, weight="bold")
    plt.yticks(size=11, weight="bold")
    plt.xlabel("$log_{2}$ fold change", size=15)
    plt.ylabel("-$log_{10}$ FDR", size=15)

    plt.legend(loc=1, bbox_to_anchor=legend_pos, frameon=False, prop={"weight": "bold"})

    if save == True:
        files = os.listdir()
        for x in range(100):
            file_pref = "volcano_" + "%02d" % (x,)
            if len([x for x in files if x.startswith(file_pref)]) == 0:
                plt.savefig(file_pref + ".png", dpi=300, bbox_inches="tight")
                plt.savefig(file_pref + ".svg", bbox_inches="tight")
                break
    elif isinstance(save, str):
        plt.savefig(save + ".png", dpi=300, bbox_inches="tight")
        plt.savefig(save + ".svg", bbox_inches="tight")

    # plt.show()
    return ax
