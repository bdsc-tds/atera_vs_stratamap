#!/usr/bin/env Rscript
# Per-cell nCount/nFeature/ratio for every gene set that meets --min-genes for this technology, +
# first_type from RCTD. Written as one small parquet per sample so the plotting scripts never need
# to reload a full Seurat object.
#
# Usage: Rscript extract_count_distributions.R --technology atera --out-dir <dir>
# (--technology must be one of CELL_LEVEL_SAMPLES$technology in plot_common.R)

suppressPackageStartupMessages({
  library(Seurat)
  library(arrow)
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--technology", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--min-genes", type = "integer", default = 200)
)))
MIN_GENES <- opt$`min-genes`
stopifnot("--technology must be one of CELL_LEVEL_SAMPLES$technology" = opt$technology %in% CELL_LEVEL_SAMPLES$technology)
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

is_stratamap <- grepl("^stratamap_", opt$technology)
cat("Technology:", opt$technology, "\n")
obj <- readRDS(qc_rds_path(opt$technology))
assay_name <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
all_genes <- rownames(obj)
cat("QCed object:", ncol(obj), "cells x", length(all_genes), "genes\n")

rctd_dir <- rctd_dir_path(opt$technology)
first_type <- rep(NA_character_, ncol(obj))
if (file.exists(file.path(rctd_dir, "spot_results.parquet"))) {
  spot_results <- read_parquet(file.path(rctd_dir, "spot_results.parquet"))
  ft <- setNames(spot_results$first_type_name, gsub("_", "-", spot_results$cell_id))
  first_type <- unname(ft[colnames(obj)])
  cat("Attached first_type:", sum(!is.na(first_type)), "/", ncol(obj), "cells\n")
}

# Trimmed to just what figures 4 (celltype_threshold_sensitivity.R, gene_set=="all" only) and 5
# (ngenes_ncells_barplot.R, "all" for the Self bar + shared_panel_genes/panel_genes for the Shared
# bar) actually read - atera_genes/xenium_5k_genes/visiumhd_genes/informative were computed by the
# original internal pipeline for other distribution plots not in this figure set.
requested_sets <- c("all")
# shared_panel_genes: the technologies it was built from (atera/visiumhd/xenium_breast100/
# stratamap); xenium_5k/xenium_biomarkers instead need "panel_genes" as their shared-panel proxy
# (ngenes_ncells_barplot.R's SHARED_PANEL_GENE_SET) since their own panels barely overlap the
# breast+100 gene list to begin with.
if (is_stratamap || opt$technology %in% c("atera", "visiumhd", "xenium_breast100")) {
  requested_sets <- c(requested_sets, "shared_panel_genes")
} else {
  requested_sets <- c(requested_sets, "panel_genes")
}

rows <- list()
for (set_name in requested_sets) {
  genes <- gene_list_for_set(set_name, all_genes)
  cat(set_name, ":", length(genes), "genes")
  if (length(genes) < MIN_GENES) { cat(" - SKIPPING (< --min-genes =", MIN_GENES, ")\n"); next }
  cat("\n")
  sub_counts <- counts[genes, , drop = FALSE]
  ncount <- Matrix::colSums(sub_counts)
  nfeature <- Matrix::colSums(sub_counts > 0)
  rows[[set_name]] <- data.frame(
    technology = opt$technology, gene_set = set_name,
    cell_id = colnames(obj), first_type = first_type,
    ncount = ncount, nfeature = nfeature, ratio = ncount / pmax(nfeature, 1),
    n_genes_in_set = length(genes), nfeature_frac = nfeature / length(genes)
  )
}

out <- do.call(rbind, rows)
out_path <- file.path(opt$`out-dir`, paste0(opt$technology, ".parquet"))
write_parquet(out, out_path)
cat("\nSaved:", out_path, "(", nrow(out), "rows)\n")
