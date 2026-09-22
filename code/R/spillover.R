#!/usr/bin/env Rscript
# Builds the spatial k-NN network and per-cell spillover metadata that figures 6/7 read
# (w2_boxplot.R, spillover_metrics.R both work off --cache-dir's cell_metadata_with_spillover.parquet
# directly, keyed by technology through the directory path - not through anything in this script's
# own output). Does a cell's RCTD-estimated secondary cell-type weight (weight_second_type) track how
# much of that type is physically nearby (neighborhood_weights_second_type, from spatial k-NN)? A
# strong relationship means the "secondary" signal looks like segmentation/diffusion spillover, not
# real co-expression or doublets. Requires the SPLIT package (https://github.com/bdsc-tds/SPLIT).
#
# Not applicable to the Chromium reference (dissociated scRNA-seq, no spatial coordinates) - this
# script is only run for technologies that have x_centroid/y_centroid.
#
# Usage: Rscript spillover.R --technology atera --qc-rds <path>_qc_seurat.rds \
#   --rctd-dir <spot_results.parquet dir> --cache-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(arrow)
  library(optparse)
})

source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/reconstruct_rctd_from_parquet.R"))

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--qc-rds", type = "character"),
  make_option("--rctd-dir", type = "character"),
  make_option("--cache-dir", type = "character"),
  make_option("--rad-pruning", type = "double", default = 15),
  make_option("--k-knn", type = "integer", default = 20)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`cache-dir`, recursive = TRUE, showWarnings = FALSE)

cat("Technology:", opt$technology, "\n")
seurat_qc <- readRDS(opt$`qc-rds`)
assay_name <- DefaultAssay(seurat_qc)

rctd <- reconstruct_rctd_from_parquet(opt$`rctd-dir`)
results_df <- rctd@results$results_df
seurat_qc <- AddMetaData(seurat_qc, results_df)
cat("spot_class table:\n"); print(table(seurat_qc$spot_class, useNA = "ifany"))

sp_coords <- as.matrix(seurat_qc@meta.data %>% select(x_centroid, y_centroid))
colnames(sp_coords) <- c("ST_1", "ST_2")
seurat_qc[["spatial"]] <- CreateDimReducObject(sp_coords, assay = assay_name, key = "ST_")

sp_nw_path <- file.path(opt$`cache-dir`, "sp_nw.rds")
if (!file.exists(sp_nw_path)) {
  sp_nw <- SPLIT::build_spatial_network(
    seurat_qc, reduction = "spatial", dims = 1:2,
    DO_prune = TRUE, rad_pruning = opt$`rad-pruning`, k_knn = opt$`k-knn`
  )
  saveRDS(sp_nw, sp_nw_path)
} else {
  sp_nw <- readRDS(sp_nw_path)
}
sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = rctd)
sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)
seurat_qc <- AddMetaData(seurat_qc, sp_neigh_df)
write_parquet(seurat_qc@meta.data, file.path(opt$`cache-dir`, "cell_metadata_with_spillover.parquet"))
cat("\nDone:", opt$`cache-dir`, "\n")
