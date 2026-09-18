#!/usr/bin/env Rscript
# Per-gene sensitivity, two normalizations - transcripts/cell/gene (mean raw count per gene across
# all QC'd cells) and transcripts/mm^2/gene (total raw count per gene divided by the sample's
# imaged tissue area). Tissue area is estimated from cell centroids (extract_spatial_metadata.R's
# output), NOT a single convex hull over all cells - a single hull would bridge any gap between
# spatially disjoint tissue pieces on the same slide, inflating area. Cell centroids are clustered
# with dbscan (eps=300um, minPts=5) first, the convex hull of each cluster taken separately, and
# their areas summed; dbscan noise points (isolated stray cells) are excluded from the area
# calculation (not from n_cells/the count matrix).
#
# Cell level only, each technology's own full POST-QC gene panel.
#
# Usage: Rscript gene_sensitivity.R --technology atera --out-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dbscan)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--technology", type = "character"),
  make_option("--out-dir", type = "character")
)))
stopifnot("--technology must be one of CELL_LEVEL_SAMPLES$technology" = opt$technology %in% CELL_LEVEL_SAMPLES$technology)
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

# Shoelace polygon area (um^2) for points already ordered around the hull (chull()'s return order).
polygon_area <- function(x, y) {
  n <- length(x)
  abs(sum(x * y[c(2:n, 1)] - x[c(2:n, 1)] * y)) / 2
}

cat("Technology:", opt$technology, "\n")
obj <- readRDS(qc_rds_path(opt$technology))
assay_name <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
n_cells <- ncol(counts)
cat("Counts:", nrow(counts), "genes x", n_cells, "cells\n")

total_transcripts <- Matrix::rowSums(counts)

spatial_f <- spatial_metadata_file_for_tech(opt$technology)
stopifnot("Spatial metadata CSV not found - run extract_spatial_metadata.R first" = file.exists(spatial_f))
sp <- read.csv(spatial_f)

db <- dbscan(as.matrix(sp[, c("x_centroid", "y_centroid")]), eps = 300, minPts = 5)
cluster_sizes <- table(db$cluster[db$cluster != 0])
cat("dbscan clusters:", length(cluster_sizes), "| sizes:", paste(sort(cluster_sizes, decreasing = TRUE), collapse = ","),
    "| noise points:", sum(db$cluster == 0), "\n")
area_um2 <- 0
for (cl in as.integer(names(cluster_sizes))) {
  pts <- sp[db$cluster == cl, c("x_centroid", "y_centroid")]
  if (nrow(pts) < 3) next
  hull_idx <- chull(pts$x_centroid, pts$y_centroid)
  area_um2 <- area_um2 + polygon_area(pts$x_centroid[hull_idx], pts$y_centroid[hull_idx])
}
area_mm2 <- area_um2 / 1e6
cat("Multi-cluster hull tissue area:", round(area_mm2, 2), "mm^2 (from", nrow(sp), "cell centroids)\n")

df <- data.frame(
  gene = names(total_transcripts),
  total_transcripts = as.numeric(total_transcripts),
  n_cells = n_cells,
  tissue_area_mm2 = area_mm2,
  transcripts_per_cell = as.numeric(total_transcripts) / n_cells,
  transcripts_per_mm2 = as.numeric(total_transcripts) / area_mm2
)
df$technology <- opt$technology
df$label <- label_for_tech(opt$technology)

out_csv <- file.path(opt$`out-dir`, paste0(opt$technology, "_gene_sensitivity.csv"))
write.csv(df, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")
