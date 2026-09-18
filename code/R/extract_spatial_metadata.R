#!/usr/bin/env Rscript
# Extracts (x_centroid, y_centroid, first_type) + technology/level/sample labels from a QCed Seurat
# object into a small CSV, so the plotting scripts never need to load a full expression matrix.
#
# Usage: Rscript extract_spatial_metadata.R --technology atera --sample breast_cancer \
#   --level cell --qc-rds <path>_qc_seurat.rds --rctd-dir <spot_results.parquet dir> --out-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
  library(optparse)
})

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character", default = NA_character_),
  make_option("--level", type = "character"),
  make_option("--qc-rds", type = "character"),
  make_option("--rctd-dir", type = "character", default = NA_character_),
  make_option("--out-dir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

cat("Technology:", opt$technology, "| Level:", opt$level, "\n")
obj <- readRDS(opt$`qc-rds`)
stopifnot("x_centroid"%in%colnames(obj@meta.data), "y_centroid"%in%colnames(obj@meta.data))

df <- data.frame(cell_id = colnames(obj), x_centroid = obj$x_centroid, y_centroid = obj$y_centroid)
if (!is.na(opt$`rctd-dir`) && file.exists(file.path(opt$`rctd-dir`, "spot_results.parquet"))) {
  spot_results <- read_parquet(file.path(opt$`rctd-dir`, "spot_results.parquet"))
  ft <- setNames(spot_results$first_type_name, gsub("_", "-", spot_results$cell_id))
  df$first_type <- unname(ft[df$cell_id])
  st <- setNames(spot_results$second_type_name, gsub("_", "-", spot_results$cell_id))
  df$second_type <- unname(st[df$cell_id])
}
df$technology <- opt$technology
df$sample <- opt$sample
df$level <- opt$level

out_path <- file.path(opt$`out-dir`, paste0("spatial_", opt$technology, "_", opt$level,
                                             ifelse(is.na(opt$sample), "", paste0("_", opt$sample)), ".csv"))
write.csv(df, out_path, row.names = FALSE)
cat("Saved", nrow(df), "rows to", out_path, "\n")
