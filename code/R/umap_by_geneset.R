#!/usr/bin/env Rscript
# Figure 2 (per-sample panels): one UMAP per cell-level sample, theme_void + aspect.ratio=1,
# colored by RCTD first_type. Each technology's own full/native panel (reductions_all.rds for
# every technology, including StrataMap - kept uniform with the ARI figure's SELF_GENE_SET rather
# than restricting StrataMap to the protein-coding+lncRNA "informative" subset).
#
# Usage: Rscript umap_by_geneset.R --out-dir <dir>

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(make_option("--out-dir", type = "character"))))
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

pal <- get_celltype_palette()

plot_bundle <- function(bundle_path, label) {
  b <- readRDS(bundle_path)
  df <- as.data.frame(b$umap)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$first_type <- factor(as.character(b$metadata$first_type), levels = names(pal))
  df <- df[!is.na(df$first_type), ]
  ggplot(df, aes(UMAP_1, UMAP_2, color = first_type)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_manual(values = pal, drop = FALSE) +
    labs(title = paste0(label, "\nN=", b$n_genes, " genes")) +
    theme_post_umap()
}

for (i in seq_len(nrow(CELL_LEVEL_SAMPLES))) {
  tech <- CELL_LEVEL_SAMPLES$technology[i]
  label <- CELL_LEVEL_SAMPLES$label[i]
  bundle_path <- reductions_path(tech, "all")
  if (!file.exists(bundle_path)) {
    cat("SKIPPING tech=", tech, "- no", bundle_path, "\n")
    next
  }
  p <- plot_bundle(bundle_path, label)
  save_post(file.path(opt$`out-dir`, paste0("umap_self_", tech, ".png")), p, "umap")
}
cat("\nDone.\n")
