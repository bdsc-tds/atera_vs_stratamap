#!/usr/bin/env Rscript
# Cell-type separation metric per technology x gene set. Clusters on the gene-set's full PCA
# embedding (k-means, k = number of distinct RCTD first_type values present) and computes ARI
# (mclust::adjustedRandIndex) against RCTD first_type. Silhouette width deliberately skipped (naive
# implementations are O(n^2), too slow at cell counts here, e.g. ~700k for Xenium 5k).
#
# Usage: Rscript separation_metric.R --technology atera \
#   --embeddings-dir <dir with reductions_*.rds> --out-csv <path>

suppressPackageStartupMessages({
  library(mclust)
  library(optparse)
})

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--embeddings-dir", type = "character"),
  make_option("--out-csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

reduction_files <- list.files(opt$`embeddings-dir`, pattern = "^reductions_.*\\.rds$", full.names = TRUE)
cat("Found", length(reduction_files), "reduction bundles in", opt$`embeddings-dir`, "\n")

results <- list()
for (f in reduction_files) {
  bundle <- readRDS(f)
  gene_set <- bundle$gene_set
  cat("\n=== Gene set:", gene_set, "===\n")

  first_type <- bundle$metadata$first_type
  if (is.null(first_type) || all(is.na(first_type))) {
    cat("SKIPPING - no RCTD first_type labels in this bundle (run after RCTD completes)\n")
    next
  }

  keep <- !is.na(first_type)
  pca <- bundle$pca[keep, , drop = FALSE]
  labels <- droplevels(as.factor(first_type[keep]))
  k <- length(levels(labels))
  cat("n_cells (annotated):", nrow(pca), "| n_types:", k, "\n")

  if (k < 2 || nrow(pca) < k * 5) {
    cat("SKIPPING - too few annotated cells/types for a meaningful clustering\n")
    next
  }

  set.seed(1)
  km <- kmeans(pca, centers = k, nstart = 10, iter.max = 100)
  ari <- adjustedRandIndex(km$cluster, labels)
  cat("ARI (kmeans on PCA vs RCTD first_type):", round(ari, 4), "\n")

  results[[gene_set]] <- data.frame(
    technology = opt$technology,
    sample = bundle$sample,
    gene_set = gene_set,
    n_genes = bundle$n_genes,
    n_cells_annotated = nrow(pca),
    n_types = k,
    ari = ari
  )
}

if (length(results) == 0) {
  cat("\nNo gene sets produced a result - nothing to save.\n")
} else {
  out_df <- do.call(rbind, results)
  dir.create(dirname(opt$`out-csv`), recursive = TRUE, showWarnings = FALSE)
  write.csv(out_df, opt$`out-csv`, row.names = FALSE)
  cat("\nSaved:", opt$`out-csv`, "\n")
  print(out_df)
}
