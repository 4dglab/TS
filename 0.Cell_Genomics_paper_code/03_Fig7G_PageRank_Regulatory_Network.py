#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_Fig7G_PageRank_Regulatory_Network.py

Purpose
-------
Calculate weighted PageRank centrality for the combined
cross-chromosomal TF-target regulatory network and generate Fig. 7G.

Input
-----
1. X-autosome-cor_sorted_v2.csv
   X-linked PAR1/escape TF -> autosomal target

2. autosome-X-cor_sorted_v2.csv
   Autosomal TF -> X-linked PAR1/escape target

Network construction
--------------------
The two regulatory directions are combined into a single directed
TF-target network:

    TF -> Target

Edge weights are defined as the absolute values of the corresponding
Spearman correlation coefficients:

    weight = abs(Correlation)

PageRank
--------
Weighted PageRank centrality is calculated for all nodes in the
combined directed network using a damping factor of 0.85.

The 20 nodes with the highest PageRank scores are shown in Fig. 7G.

Outputs
-------
1. pagerank_scores.csv
   PageRank scores for all nodes.

2. pagerank_scores_top20.csv
   Top 20 PageRank-ranked nodes.

3. Fig7G_PageRank_Top20_Nodes.pdf
   Fig. 7G.
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx


# ============================================================
# 1. Paths
# ============================================================

DATA_DIR = Path("Source data")
OUTPUT_DIR = Path("Figures")
RESULT_DIR = Path("Results")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)

X_TO_AUTOSOME_FILE = (
    DATA_DIR /
    "X-autosome-cor_sorted_v2.csv"
)

AUTOSOME_TO_X_FILE = (
    DATA_DIR /
    "autosome-X-cor_sorted_v2.csv"
)


# ============================================================
# 2. Load TF-target interactions
# ============================================================

# X-linked PAR1/escape TF -> autosomal target
x_to_autosome = pd.read_csv(
    X_TO_AUTOSOME_FILE
)

# Autosomal TF -> X-linked PAR1/escape target
autosome_to_x = pd.read_csv(
    AUTOSOME_TO_X_FILE
)


# ============================================================
# 3. Combine the two regulatory directions
# ============================================================

x_to_autosome["Regulation"] = "X-to-Autosome"

autosome_to_x["Regulation"] = "Autosome-to-X"


combined_data = pd.concat(
    [
        x_to_autosome,
        autosome_to_x
    ],
    ignore_index=True
)


# Ensure correlation values are numeric
combined_data["Correlation"] = pd.to_numeric(
    combined_data["Correlation"],
    errors="coerce"
)

combined_data = combined_data.dropna(
    subset=[
        "TF",
        "Target",
        "Correlation"
    ]
).copy()


# ============================================================
# 4. Construct directed weighted regulatory network
# ============================================================

G = nx.DiGraph()


for _, row in combined_data.iterrows():

    tf = row["TF"]
    target = row["Target"]

    # Regulatory direction:
    #
    #       TF  ->  Target
    #
    # Edge strength is represented by the absolute
    # Spearman correlation coefficient.

    weight = abs(
        row["Correlation"]
    )

    G.add_edge(
        tf,
        target,
        weight=weight
    )


# ============================================================
# 5. Calculate weighted PageRank
# ============================================================

pagerank_scores = nx.pagerank(
    G,
    alpha=0.85,
    weight="weight"
)


pagerank_df = pd.DataFrame(
    pagerank_scores.items(),
    columns=[
        "Node",
        "PageRank"
    ]
)


pagerank_df = (
    pagerank_df
    .sort_values(
        by="PageRank",
        ascending=False
    )
    .reset_index(drop=True)
)


# ============================================================
# 6. Save PageRank scores for all network nodes
# ============================================================

pagerank_df.to_csv(
    RESULT_DIR /
    "pagerank_scores.csv",
    index=False
)


# ============================================================
# 7. Select Top 20 PageRank nodes
# ============================================================

pagerank_top20 = (
    pagerank_df
    .head(20)
    .copy()
)


pagerank_top20.to_csv(
    RESULT_DIR /
    "pagerank_scores_top20.csv",
    index=False
)


# ============================================================
# 8. Generate Fig. 7G
# ============================================================

plt.figure(
    figsize=(8, 6)
)


sns.barplot(
    data=pagerank_top20,
    x="PageRank",
    y="Node",
    hue="Node",
    palette="Reds_r",
    legend=False
)


plt.xlabel(
    "PageRank score"
)

plt.ylabel(
    "Node"
)

plt.title(
    "Top 20 Nodes in the Regulatory Network"
)


plt.savefig(
    OUTPUT_DIR /
    "Fig7G_PageRank_Top20_Nodes.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 9. Report network statistics
# ============================================================

print("\nPageRank analysis summary")
print("=" * 60)

print(
    f"X-to-autosome interactions: "
    f"{len(x_to_autosome)}"
)

print(
    f"Autosome-to-X interactions: "
    f"{len(autosome_to_x)}"
)

print(
    f"Combined input interactions: "
    f"{len(combined_data)}"
)

print(
    f"Network nodes: "
    f"{G.number_of_nodes()}"
)

print(
    f"Network edges: "
    f"{G.number_of_edges()}"
)


print("\nTop 20 PageRank nodes:")
print(
    pagerank_top20.to_string(
        index=False
    )
)


print("\nOutputs:")
print(
    RESULT_DIR /
    "pagerank_scores.csv"
)

print(
    RESULT_DIR /
    "pagerank_scores_top20.csv"
)

print(
    OUTPUT_DIR /
    "Fig7G_PageRank_Top20_Nodes.pdf"
)