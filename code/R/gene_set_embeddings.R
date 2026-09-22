#!/usr/bin/env Rscript
# Per-technology, per-gene-set PCA/UMAP. HVGs are computed WITHIN each gene-set subset (subset the
# Seurat object first, THEN FindVariableFeatures) - not genome-wide HVGs intersected afterward,
# which silently distorts the embedding for small gene sets.
#
# Gene sets (see plot_common.R's gene_list_for_set(); which ones apply to a given technology is
# decided by the CALLING script, not hardcoded here):
#   all                - every gene in the QCed object
#   panel_genes        - intersect with reference/breast100_gene_list.csv (breast+100 addon panel)
#   shared_panel_genes - intersect with reference/shared_panel_gene_list.csv (single literal gene
#                        list shared across every whole-transcriptome technology's QCed data)
# A gene set is skipped (with a loud log line, not silently) if the intersection has too few genes
# to be meaningful (--min-genes, default 200).
#
# --min-count-in-subset / --npcs / --label-suffix: a small gene set (e.g. panel_genes) restricted
# from a whole-transcriptome technology's full panel can leave many cells near-empty even though
# they passed the object's overall QC threshold - the symptom is a starburst/spoke UMAP artifact.
# Rerunning that gene set with --min-count-in-subset (drop near-empty cells WITHIN the subset) and a
# lower --npcs (fewer PCs for a smaller/sparser gene set) produces a cleaner embedding;
# --label-suffix (e.g. "_filtered") saves it alongside the original instead of overwriting it.
#
# Usage: Rscript gene_set_embeddings.R --technology atera --sample breast_cancer \
#   --qc-rds <path>_qc_seurat.rds --gene-sets all,panel_genes \
#   --out-dir <dir> [--rctd-dir <spot_results.parquet dir>]
#   [--min-count-in-subset 20 --npcs 15 --label-suffix _filtered] [--umap-init random]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(arrow)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character"),
  make_option("--qc-rds", type = "character"),
  make_option("--gene-sets", type = "character", default = "all,panel_genes"),
  make_option("--rctd-dir", type = "character", default = NA),
  make_option("--out-dir", type = "character"),
  make_option("--min-count-in-subset", type = "integer", default = 0),
  make_option("--npcs", type = "integer", default = 50),
  make_option("--label-suffix", type = "character", default = ""),
  make_option("--umap-init", type = "character", default = "spectral",
              help = "'spectral' (Seurat::RunUMAP default) or 'random' - use 'random' when spectral init segfaults RSpectra at very large cell counts (~500k+)"),
  make_option("--min-genes", type = "integer", default = 200)
)
opt <- parse_args(OptionParser(option_list = option_list))
requested_sets <- strsplit(opt$`gene-sets`, ",")[[1]]
MIN_GENES <- opt$`min-genes`

cat("Technology:", opt$technology, "| Sample:", opt$sample, "\n")
obj <- readRDS(opt$`qc-rds`)
assay_name <- DefaultAssay(obj)
all_genes <- rownames(obj)
cat("QCed object:", ncol(obj), "cells x", length(all_genes), "genes\n")

if (!is.na(opt$`rctd-dir`) && file.exists(file.path(opt$`rctd-dir`, "spot_results.parquet"))) {
  spot_results <- read_parquet(file.path(opt$`rctd-dir`, "spot_results.parquet"))
  spot_results$cell_id <- gsub("_", "-", spot_results$cell_id)
  rownames(spot_results) <- spot_results$cell_id
  meta_cols <- c("spot_class", "first_type_name", "second_type_name")
  spot_results <- spot_results[colnames(obj), meta_cols]
  colnames(spot_results) <- c("spot_class", "first_type", "second_type")
  obj <- AddMetaData(obj, spot_results)
  cat("Attached RCTD metadata:", sum(!is.na(obj$first_type)), "/", ncol(obj), "cells annotated\n")
} else {
  cat("No RCTD dir supplied or spot_results.parquet not found - proceeding without cell-type metadata\n")
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

for (set_name in requested_sets) {
  cat("\n=== Gene set:", set_name, "===\n")
  genes <- gene_list_for_set(set_name, all_genes)
  cat("Genes in set (intersected with object):", length(genes), "\n")
  if (length(genes) < MIN_GENES) {
    cat("SKIPPING", set_name, "- only", length(genes), "genes, below --min-genes =", MIN_GENES, "\n")
    next
  }

  sub_obj <- subset(obj, features = genes)

  if (opt$`min-count-in-subset` > 0) {
    subset_counts <- Matrix::colSums(GetAssayData(sub_obj, assay = assay_name, layer = "counts"))
    n_before <- ncol(sub_obj)
    sub_obj <- subset(sub_obj, cells = colnames(sub_obj)[subset_counts >= opt$`min-count-in-subset`])
    cat("Filtered cells with <", opt$`min-count-in-subset`, "counts within this", length(genes),
        "-gene subset:", n_before, "->", ncol(sub_obj), "cells\n")
  }

  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  nfeatures <- min(2000, nrow(sub_obj) - 1)
  sub_obj <- FindVariableFeatures(sub_obj, nfeatures = nfeatures, verbose = FALSE)
  cat("HVGs within", set_name, "set:", length(VariableFeatures(sub_obj)), "\n")
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  npcs <- min(opt$npcs, ncol(sub_obj) - 1, length(VariableFeatures(sub_obj)) - 1)
  sub_obj <- RunPCA(sub_obj, npcs = npcs, verbose = FALSE)

  # Seurat::RunUMAP hardcodes uwot's spectral graph-Laplacian init internally - at very large cell
  # counts this can segfault inside RSpectra::eigs_sym (not an OOM). --umap-init random bypasses
  # Seurat::RunUMAP entirely and calls uwot::umap() directly on the PCA embedding with random init
  # instead, same neighbor/metric/min_dist settings as Seurat's own defaults.
  if (opt$`umap-init` == "random") {
    cat("Using uwot::umap() directly with init='random' (bypassing Seurat::RunUMAP's spectral init)\n")
    umap_coords <- uwot::umap(Embeddings(sub_obj, "pca")[, 1:npcs], n_neighbors = 30L, metric = "cosine",
                               min_dist = 0.3, spread = 1, init = "random", seed = 42L, verbose = FALSE)
    rownames(umap_coords) <- colnames(sub_obj)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  } else {
    sub_obj <- RunUMAP(sub_obj, dims = 1:npcs, verbose = FALSE)
    umap_coords <- Embeddings(sub_obj, "umap")
  }

  out_set_name <- paste0(set_name, opt$`label-suffix`)
  meta_cols_present <- intersect(c("first_type", "second_type", "spot_class", paste0("nCount_", assay_name), paste0("nFeature_", assay_name)), colnames(sub_obj@meta.data))
  bundle <- list(
    pca = Embeddings(sub_obj, "pca"),
    umap = umap_coords,
    metadata = sub_obj@meta.data[, meta_cols_present, drop = FALSE],
    hvg = VariableFeatures(sub_obj),
    gene_set = out_set_name,
    n_genes = length(genes),
    technology = opt$technology,
    sample = opt$sample
  )
  out_path <- file.path(opt$`out-dir`, paste0("reductions_", out_set_name, ".rds"))
  saveRDS(bundle, out_path)
  cat("Saved:", out_path, "\n")
  rm(sub_obj); gc()
}

cat("\nDone.\n")
