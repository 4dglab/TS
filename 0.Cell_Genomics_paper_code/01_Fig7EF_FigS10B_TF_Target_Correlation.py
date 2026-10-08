#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_Fig7EF_FigS10B_TF_Target_Correlation.py

Purpose
-------
Generate the TF-target correlation plots shown in:
    - Fig. 7E-F
    - Fig. S10B

Input
-----
1. X-autosome-cor_sorted_v2.csv
   X-linked TF -> autosomal target

2. autosome-X-cor_sorted_v2.csv
   Autosomal TF -> X-linked PAR1/escape target

Analysis
--------
1. Combine the two regulatory directions.
2. Rank all TF-target interactions globally by the absolute value of
   the Spearman correlation coefficient.
3. Annotate the five strongest positive and five strongest negative
   correlations within each regulatory direction.
4. Fig. 7E-F:
   visualize each regulatory direction separately and perform a
   two-sided Wilcoxon signed-rank test against zero.
5. Fig. S10B:
   visualize the two regulatory directions jointly, together with
   marginal correlation distributions.

Outputs
-------
Fig7EF_TF_Target_Correlation_by_Direction.pdf
FigS10B_Global_TF_Target_Correlation.pdf
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
from adjustText import adjust_text


# ============================================================
# 1. Paths
# ============================================================

DATA_DIR = Path("Source data")
OUTPUT_DIR = Path("Figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

X_TO_AUTOSOME_FILE = DATA_DIR / "X-autosome-cor_sorted_v2.csv"
AUTOSOME_TO_X_FILE = DATA_DIR / "autosome-X-cor_sorted_v2.csv"


# ============================================================
# 2. Plot settings
# ============================================================

sns.set_theme(
    style="whitegrid",
    font_scale=1.2
)

plt.rcParams["font.family"] = "Arial"

COLOR_DICT = {
    "Autosome-to-X": "#6B5B95",
    "X-to-Autosome": "#FF6F61"
}

# Explicit order used in Fig. 7E-F
TYPE_ORDER = [
    "Autosome-to-X",
    "X-to-Autosome"
]


# ============================================================
# 3. Load TF-target correlation data
# ============================================================

# X-linked TF -> autosomal target
x_to_autosome = pd.read_csv(
    X_TO_AUTOSOME_FILE
)

x_to_autosome["type"] = "X-to-Autosome"


# Autosomal TF -> X-linked PAR1/escape target
autosome_to_x = pd.read_csv(
    AUTOSOME_TO_X_FILE
)

autosome_to_x["type"] = "Autosome-to-X"


# Combine both regulatory directions
combined_df = pd.concat(
    [x_to_autosome, autosome_to_x],
    ignore_index=True
)


# ============================================================
# 4. Prepare correlation data
# ============================================================

combined_df["Correlation"] = pd.to_numeric(
    combined_df["Correlation"],
    errors="coerce"
)

combined_df = combined_df.dropna(
    subset=["Correlation"]
).copy()


# Rank all interactions globally by absolute Spearman correlation
combined_df["Abs_Correlation"] = combined_df["Correlation"].abs()

combined_df = (
    combined_df
    .sort_values(
        "Abs_Correlation",
        ascending=False
    )
    .reset_index(drop=True)
)

# This corresponds to the "index" variable in the original analysis
combined_df["Rank"] = np.arange(len(combined_df))


# ============================================================
# 5. Select interactions for annotation
# ============================================================

# For each regulatory direction:
#   - five strongest positive correlations
#   - five strongest negative correlations

highlight_list = []

for regulation_type in TYPE_ORDER:

    subset = combined_df[
        combined_df["type"] == regulation_type
    ]

    top5_positive = subset.nlargest(
        5,
        "Correlation"
    )

    top5_negative = subset.nsmallest(
        5,
        "Correlation"
    )

    highlight_list.extend(
        [top5_positive, top5_negative]
    )


highlight_points = pd.concat(
    highlight_list,
    ignore_index=True
)


# ============================================================
# 6. Fig. S10B
#    Global distribution of TF-target correlations
# ============================================================

g = sns.jointplot(
    data=combined_df,
    x="Rank",
    y="Correlation",
    hue="type",
    hue_order=TYPE_ORDER,
    kind="scatter",
    height=8,
    alpha=0.6,
    palette=COLOR_DICT
)

ax = g.ax_joint


# Annotate top/bottom five interactions from each direction
texts = []

for _, row in highlight_points.iterrows():

    label = f"{row['TF']}-{row['Target']}"

    texts.append(
        ax.text(
            row["Rank"],
            row["Correlation"],
            label,
            fontsize=10,
            weight="bold",
            color=COLOR_DICT[row["type"]]
        )
    )


adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(
        arrowstyle="->",
        color="gray",
        lw=0.5
    )
)


# Labels
g.ax_joint.set_xlabel("Data points")
g.ax_joint.set_ylabel("Spearman correlation")

plt.suptitle(
    "Global Distribution of TF-Target Correlations",
    y=1.02
)


# Save Fig. S10B
plt.savefig(
    OUTPUT_DIR / "FigS10B_Global_TF_Target_Correlation.pdf",
    dpi=1200,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 7. Fig. 7E-F
#    Correlation distributions by regulatory direction
# ============================================================

g = sns.FacetGrid(
    combined_df,
    col="type",
    col_order=TYPE_ORDER,
    height=6,
    aspect=1.2,
    sharey=False
)


# ------------------------------------------------------------
# Scatter plot
# ------------------------------------------------------------

def scatter_with_color(data, **kwargs):

    regulation_type = data["type"].iloc[0]

    sns.scatterplot(
        data=data,
        x="Rank",
        y="Correlation",
        alpha=0.7,
        edgecolor="white",
        linewidth=0.5,
        s=80,
        color=COLOR_DICT[regulation_type]
    )


g.map_dataframe(
    scatter_with_color
)


# ------------------------------------------------------------
# 2D KDE contours
# ------------------------------------------------------------

g.map_dataframe(
    sns.kdeplot,
    x="Rank",
    y="Correlation",
    levels=5,
    color="black",
    alpha=0.5
)


# ------------------------------------------------------------
# Reference line at correlation = 0
# ------------------------------------------------------------

for ax in g.axes.flat:

    ax.axhline(
        y=0,
        linestyle="--",
        color="red",
        alpha=0.5
    )


# ------------------------------------------------------------
# Median and Wilcoxon signed-rank test
# ------------------------------------------------------------

def annotate_stats(data, **kwargs):

    correlations = (
        data["Correlation"]
        .dropna()
        .to_numpy()
    )

    median = np.median(correlations)

    if len(correlations) > 1:

        _, p_value = stats.wilcoxon(
            correlations,
            alternative="two-sided"
        )

    else:

        p_value = np.nan


    ax = plt.gca()

    ax.text(
        0.70,
        0.90,
        f"Median = {median:.2f}\n"
        f"p = {p_value:.2e}",
        transform=ax.transAxes,
        fontsize=12
    )


g.map_dataframe(
    annotate_stats
)


# ------------------------------------------------------------
# Axis ranges
# ------------------------------------------------------------

for ax in g.axes.flat:

    ax.set_ylim(
        -0.10,
        0.11
    )

    ax.set_xlim(
        -190,
        1550
    )


# ------------------------------------------------------------
# Annotate top/bottom interactions
# ------------------------------------------------------------

for ax, regulation_type in zip(
    g.axes.flat,
    TYPE_ORDER
):

    subset = highlight_points[
        highlight_points["type"] == regulation_type
    ]

    texts = []

    for _, row in subset.iterrows():

        label = (
            f"{row['TF']}-"
            f"{row['Target']}"
        )

        texts.append(
            ax.text(
                row["Rank"],
                row["Correlation"],
                label,
                fontsize=9,
                weight="bold",
                color=COLOR_DICT[regulation_type]
            )
        )


    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(
            arrowstyle="->",
            color="gray",
            lw=0.5
        ),
        only_move={
            "points": "y",
            "text": "xy"
        }
    )


# ------------------------------------------------------------
# Panel titles
# ------------------------------------------------------------

g.axes.flat[0].set_title(
    "TF on autosome and target on X"
)

g.axes.flat[1].set_title(
    "TF on X and target on autosome"
)


# ------------------------------------------------------------
# Axis labels
# ------------------------------------------------------------

g.set_axis_labels(
    "Data points",
    "Spearman correlation"
)


# Add y-axis label to the right panel,
# matching the original Fig. 7F layout
right_ax = g.axes.flat[1]

right_ax.annotate(
    "Spearman correlation",
    xy=(-0.17, 0.5),
    xycoords="axes fraction",
    ha="center",
    va="center",
    rotation=90,
    fontsize=14
)


plt.subplots_adjust(
    wspace=0.25
)


# ------------------------------------------------------------
# Save Fig. 7E-F
# ------------------------------------------------------------

plt.savefig(
    OUTPUT_DIR / "Fig7EF_TF_Target_Correlation_by_Direction.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 8. Report summary statistics
# ============================================================

print("\nAnalysis summary")
print("=" * 60)

print(
    f"Total TF-target interactions: "
    f"{len(combined_df)}"
)

for regulation_type in TYPE_ORDER:

    subset = combined_df[
        combined_df["type"] == regulation_type
    ]

    correlations = subset["Correlation"].dropna()

    _, p_value = stats.wilcoxon(
        correlations,
        alternative="two-sided"
    )

    print(
        f"\n{regulation_type}"
        f"\n  Interactions: {len(subset)}"
        f"\n  Median Spearman correlation: "
        f"{correlations.median():.6f}"
        f"\n  Wilcoxon signed-rank P: "
        f"{p_value:.6e}"
    )


print("\nFigures saved to:")
print(
    OUTPUT_DIR /
    "FigS10B_Global_TF_Target_Correlation.pdf"
)

print(
    OUTPUT_DIR /
    "Fig7EF_TF_Target_Correlation_by_Direction.pdf"
)