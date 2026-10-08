#!/usr/bin/env Rscript

# ============================================================
# 00_Preprocess_TF_Target_Edges.R
#
# Description:
#   Annotate cross-chromosomal TF-target interactions
#   with previously calculated correlation information.
#
# Inputs:
#   edge_info.tsv
#   autosome-regulate-X.csv
#   X-regulate-autosome.csv
#
# Outputs:
#   autosome-regulate-X-cor.csv
#   X-regulate-autosome-cor.csv
#
# Note:
#   Correlations are obtained from an existing edge table.
#   No correlation coefficients are recalculated.
# ============================================================

library(readr)
library(dplyr)


# ============================================================
# 1. Load correlation information
# ============================================================

edge.mat <- read_tsv("edge_info.tsv")


# ============================================================
# 2. Define correlation annotation function
# ============================================================

add_cor <- function(df) {

  df <- df %>%
    left_join(
      edge.mat,
      by = c(
        "TF" = "Source",
        "target" = "Target"
      )
    )

  return(df)
}


# ============================================================
# 3. Autosome TF -> X chromosome target
# ============================================================

auto_x <- read_csv("autosome-regulate-X.csv")

auto_x <- add_cor(auto_x)


# ============================================================
# 4. X chromosome TF -> autosomal target
# ============================================================

x_auto <- read_csv("X-regulate-autosome.csv")

x_auto <- add_cor(x_auto)


# ============================================================
# 5. Export annotated TF-target interactions
# ============================================================

write_csv(
  auto_x,
  "autosome-regulate-X-cor.csv"
)

write_csv(
  x_auto,
  "X-regulate-autosome-cor.csv"
)