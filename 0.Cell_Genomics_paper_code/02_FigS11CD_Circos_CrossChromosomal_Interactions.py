#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_FigS11CD_Circos_CrossChromosomal_Interactions.py

Purpose
-------
Generate Circos plots showing the genomic distribution of
cross-chromosomal TF-target regulatory interactions.

Panels
------
Fig. S11C:
    Autosomal TF -> X-linked PAR1/escape target

Fig. S11D:
    X-linked PAR1/escape TF -> autosomal target

Input
-----
1. autosome-X-cor_sorted_v2.csv
2. X-autosome-cor_sorted_v2.csv
3. gene_chromosome_map_updated.csv

Genome assembly
---------------
hg38

Visualization
-------------
Positive Spearman correlation: red
Negative Spearman correlation: blue
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from pycirclize import Circos
from pycirclize.utils import load_eukaryote_example_dataset


# ============================================================
# 1. Paths
# ============================================================

DATA_DIR = Path("Source data")
OUTPUT_DIR = Path("Figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AUTOSOME_TO_X_FILE = DATA_DIR / "autosome-X-cor_sorted_v2.csv"
X_TO_AUTOSOME_FILE = DATA_DIR / "X-autosome-cor_sorted_v2.csv"

GENE_COORD_FILE = DATA_DIR / "gene_chromosome_map_updated.csv"


# ============================================================
# 2. Load hg38 chromosome and cytoband information
# ============================================================

chr_bed_file, cytoband_file, _ = (
    load_eukaryote_example_dataset("hg38")
)


# ============================================================
# 3. Load regulatory interactions
# ============================================================

# Autosomal TF -> X-linked PAR1/escape target
autosome_to_x = pd.read_csv(
    AUTOSOME_TO_X_FILE
)

# X-linked PAR1/escape TF -> autosomal target
x_to_autosome = pd.read_csv(
    X_TO_AUTOSOME_FILE
)


# ============================================================
# 4. Load gene genomic coordinates
# ============================================================

gene_coords = pd.read_csv(
    GENE_COORD_FILE
)


# ============================================================
# 5. Standardize gene symbols
# ============================================================

gene_coords["gene"] = (
    gene_coords["gene"]
    .astype(str)
    .str.upper()
)

for df in [autosome_to_x, x_to_autosome]:

    df["TF"] = (
        df["TF"]
        .astype(str)
        .str.upper()
    )

    df["Target"] = (
        df["Target"]
        .astype(str)
        .str.upper()
    )


# ============================================================
# 6. Gene-coordinate lookup
# ============================================================

def get_gene_region(gene):
    """
    Return genomic coordinates for a gene in pyCirclize format.

    Returns
    -------
    tuple
        (chromosome, start, end)

    None
        If the gene is absent from the coordinate table.
    """

    row = gene_coords[
        gene_coords["gene"] == gene
    ]

    if row.empty:
        return None

    chrom = str(
        row.iloc[0]["chrom"]
    )

    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    start = int(
        row.iloc[0]["start"]
    )

    end = int(
        row.iloc[0]["end"]
    )

    return chrom, start, end


# ============================================================
# 7. Circos plotting function
# ============================================================

def plot_regulatory_circos(
    interactions,
    output_file
):
    """
    Plot cross-chromosomal TF-target interactions.

    Parameters
    ----------
    interactions : pandas.DataFrame
        Must contain TF, Target, and Correlation columns.

    output_file : pathlib.Path
        Output PDF path.
    """

    # --------------------------------------------------------
    # Initialize hg38 Circos plot
    # --------------------------------------------------------

    circos = Circos.initialize_from_bed(
        chr_bed_file,
        space=3
    )


    # --------------------------------------------------------
    # Add cytoband tracks
    # --------------------------------------------------------

    circos.add_cytoband_tracks(
        (95, 100),
        cytoband_file
    )


    # --------------------------------------------------------
    # Chromosome labels and genomic coordinates
    # --------------------------------------------------------

    for sector in circos.sectors:

        sector.text(
            sector.name,
            r=120,
            size=10,
            color="#333333"
        )

        sector.get_track(
            "cytoband"
        ).xticks_by_interval(
            40_000_000,
            label_size=8,
            label_orientation="vertical",
            label_formatter=lambda v:
                f"{v / 1_000_000:.0f} Mb"
        )


    # --------------------------------------------------------
    # Add TF-target links
    # --------------------------------------------------------

    n_plotted = 0
    missing_genes = set()

    for _, row in interactions.iterrows():

        tf_region = get_gene_region(
            row["TF"]
        )

        target_region = get_gene_region(
            row["Target"]
        )


        if tf_region is None:

            missing_genes.add(
                row["TF"]
            )

            continue


        if target_region is None:

            missing_genes.add(
                row["Target"]
            )

            continue


        # Positive correlation = red
        # Negative correlation = blue

        link_color = (
            "red"
            if row["Correlation"] > 0
            else "blue"
        )


        circos.link(
            tf_region,
            target_region,
            color=link_color
        )

        n_plotted += 1


    # --------------------------------------------------------
    # Render and save
    # --------------------------------------------------------

    circos.plotfig()

    plt.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()


    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    print(
        f"\nSaved: {output_file}"
    )

    print(
        f"Input interactions: "
        f"{len(interactions)}"
    )

    print(
        f"Interactions plotted: "
        f"{n_plotted}"
    )

    if missing_genes:

        print(
            f"Genes without genomic coordinates: "
            f"{len(missing_genes)}"
        )

        print(
            ", ".join(
                sorted(missing_genes)
            )
        )


# ============================================================
# 8. Fig. S11C
#    Autosomal TF -> X-linked target
# ============================================================

plot_regulatory_circos(
    interactions=autosome_to_x,
    output_file=(
        OUTPUT_DIR /
        "FigS11C_Autosome_to_X_Circos.pdf"
    )
)


# ============================================================
# 9. Fig. S11D
#    X-linked TF -> autosomal target
# ============================================================

plot_regulatory_circos(
    interactions=x_to_autosome,
    output_file=(
        OUTPUT_DIR /
        "FigS11D_X_to_Autosome_Circos.pdf"
    )
)