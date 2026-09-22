#!/usr/bin/env Rscript
# Sensitivity comparison across technologies: nCount, nGenes (nFeature), and nCount/nGenes ratio,
# per gene set, across the whole dataset (feeds figure 5's --overall-csv, via sensitivity_synthesis.R).
#
# --min-count-in-subset (default 10, applied to every gene set except "all"): matches both
# gene_set_embeddings.R's own --min-count-in-subset for panel_genes/shared_panel_genes (STAGE 4b/4c)
# AND the self/"all" panel's own QC floor (qc_seurat.R's --min-counts default of 10) - the SAME
# cell-inclusion threshold everywhere, so figure 5's n_cells/ARI panels (from the embeddings
# pipeline) and its nCount/nCount-per-nGenes panels (from this script) describe the same cells for
# a given (technology, gene set), not two different populations under the same "shared panel genes"
# label. Never applied to the 'all' gene set (its embedding isn't filtered either - it's already
# exactly the QC'd object).
#
# Usage: Rscript sensitivity_by_geneset_celltype.R --technology atera \
#   --qc-rds <path>_qc_seurat.rds --gene-sets all,panel_genes --out-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character", default = NA_character_),
  make_option("--qc-rds", type = "character"),
  make_option("--gene-sets", type = "character", default = "all,panel_genes"),
  make_option("--min-count-in-subset", type = "integer", default = 10,
              help = "drop cells with fewer than this many counts WITHIN the gene-set subset before summarising - never applied to the 'all' gene set"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
requested_sets <- strsplit(opt$`gene-sets`, ",")[[1]]

cat("Technology:", opt$technology, "\n")
obj <- readRDS(opt$`qc-rds`)
assay_name <- DefaultAssay(obj)
counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
all_genes <- rownames(obj)

summarise_vec <- function(x) {
  x <- x[is.finite(x)]
  data.frame(n = length(x), mean = mean(x), median = median(x), sd = sd(x),
             p25 = quantile(x, 0.25), p75 = quantile(x, 0.75), min = min(x), max = max(x))
}

overall_rows <- list()

for (set_name in requested_sets) {
  genes <- gene_list_for_set(set_name, all_genes)
  cat(set_name, ":", length(genes), "genes\n")
  if (length(genes) == 0) { cat("SKIPPING - 0 genes\n"); next }

  sub_counts <- counts[genes, , drop = FALSE]
  nCount <- Matrix::colSums(sub_counts)
  nGenes <- Matrix::colSums(sub_counts > 0)

  if (set_name != "all" && opt$`min-count-in-subset` > 0) {
    keep <- nCount >= opt$`min-count-in-subset`
    cat("Filtered cells with <", opt$`min-count-in-subset`, "counts within this", length(genes),
        "-gene subset:", length(nCount), "->", sum(keep), "cells\n")
    nCount <- nCount[keep]; nGenes <- nGenes[keep]
  }
  ratio <- ifelse(nGenes > 0, nCount / nGenes, NA_real_)

  ov <- rbind(
    cbind(metric = "nCount", summarise_vec(nCount)),
    cbind(metric = "nGenes", summarise_vec(nGenes)),
    cbind(metric = "nCount_per_nGenes", summarise_vec(ratio))
  )
  ov$technology <- opt$technology; ov$sample <- opt$sample; ov$gene_set <- set_name; ov$n_genes_in_set <- length(genes)
  overall_rows[[set_name]] <- ov
}

overall_df <- do.call(rbind, overall_rows)
write.csv(overall_df, file.path(opt$`out-dir`, paste0("sensitivity_overall_", opt$technology, ".csv")), row.names = FALSE)
cat("\nSaved overall summary:", nrow(overall_df), "rows\n")
cat("\nDone:", opt$`out-dir`, "\n")
