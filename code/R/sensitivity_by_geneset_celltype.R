#!/usr/bin/env Rscript
# Sensitivity comparison across technologies: nCount, nGenes (nFeature), and nCount/nGenes ratio,
# per gene set, computed both across the whole dataset and per RCTD cell type.
#
# Usage: Rscript sensitivity_by_geneset_celltype.R --technology atera \
#   --qc-rds <path>_qc_seurat.rds --rctd-dir <spot_results.parquet dir> \
#   --gene-sets all,atera_genes,panel_genes --out-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character", default = NA_character_),
  make_option("--qc-rds", type = "character"),
  make_option("--rctd-dir", type = "character", default = NA_character_,
              help = "optional - if omitted, only the whole-dataset summary is produced"),
  make_option("--gene-sets", type = "character", default = "all,atera_genes,panel_genes"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
requested_sets <- strsplit(opt$`gene-sets`, ",")[[1]]

cat("Technology:", opt$technology, "\n")
obj <- readRDS(opt$`qc-rds`)
assay_name <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
all_genes <- rownames(obj)

first_type <- rep(NA_character_, ncol(obj))
names(first_type) <- colnames(obj)
if (!is.na(opt$`rctd-dir`)) {
  spot_results <- read_parquet(file.path(opt$`rctd-dir`, "spot_results.parquet"))
  ft <- setNames(spot_results$first_type_name, gsub("_", "-", spot_results$cell_id))
  matched <- intersect(colnames(obj), names(ft))
  first_type[matched] <- unname(ft[matched])
  cat(sum(!is.na(first_type)), "/", length(first_type), "cells have an RCTD cell type\n")
}

summarise_vec <- function(x) {
  x <- x[is.finite(x)]
  data.frame(n = length(x), mean = mean(x), median = median(x), sd = sd(x),
             p25 = quantile(x, 0.25), p75 = quantile(x, 0.75), min = min(x), max = max(x))
}

overall_rows <- list()
celltype_rows <- list()

for (set_name in requested_sets) {
  genes <- gene_list_for_set(set_name, all_genes)
  cat(set_name, ":", length(genes), "genes\n")
  if (length(genes) == 0) { cat("SKIPPING - 0 genes\n"); next }

  sub_counts <- counts[genes, , drop = FALSE]
  nCount <- Matrix::colSums(sub_counts)
  nGenes <- Matrix::colSums(sub_counts > 0)
  ratio <- ifelse(nGenes > 0, nCount / nGenes, NA_real_)

  ov <- rbind(
    cbind(metric = "nCount", summarise_vec(nCount)),
    cbind(metric = "nGenes", summarise_vec(nGenes)),
    cbind(metric = "nCount_per_nGenes", summarise_vec(ratio))
  )
  ov$technology <- opt$technology; ov$sample <- opt$sample; ov$gene_set <- set_name; ov$n_genes_in_set <- length(genes)
  overall_rows[[set_name]] <- ov

  if (!is.na(opt$`rctd-dir`)) {
    for (ct in na.omit(unique(first_type))) {
      idx <- which(first_type == ct)
      ct_df <- rbind(
        cbind(metric = "nCount", summarise_vec(nCount[idx])),
        cbind(metric = "nGenes", summarise_vec(nGenes[idx])),
        cbind(metric = "nCount_per_nGenes", summarise_vec(ratio[idx]))
      )
      ct_df$technology <- opt$technology; ct_df$sample <- opt$sample; ct_df$gene_set <- set_name
      ct_df$cell_type <- ct; ct_df$n_genes_in_set <- length(genes)
      celltype_rows[[paste(set_name, ct)]] <- ct_df
    }
  }
}

overall_df <- do.call(rbind, overall_rows)
write.csv(overall_df, file.path(opt$`out-dir`, paste0("sensitivity_overall_", opt$technology, ".csv")), row.names = FALSE)
cat("\nSaved overall summary:", nrow(overall_df), "rows\n")

if (length(celltype_rows) > 0) {
  celltype_df <- do.call(rbind, celltype_rows)
  write.csv(celltype_df, file.path(opt$`out-dir`, paste0("sensitivity_by_celltype_", opt$technology, ".csv")), row.names = FALSE)
  cat("Saved per-cell-type summary:", nrow(celltype_df), "rows\n")
}
cat("\nDone:", opt$`out-dir`, "\n")
