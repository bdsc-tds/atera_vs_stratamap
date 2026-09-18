#!/usr/bin/env Rscript
# Figure 2 (grid): single row-grid composite of the 8 "self" UMAPs (each technology's own
# full/native panel - same gene-set-per-technology logic as umap_by_geneset.R's GROUPS$self),
# companion to spatial_grid.R/composition_grid.R which do the same row-grid treatment for the
# spatial and composition figures. Every panel is aspect=1 (theme_post_umap() hardcodes it).
#
# Usage: Rscript umap_self_grid.R --out-dir <dir>
#   [--panel-width 4] [--panel-height 4] [--format png|pdf] [--dpi 300] [--out-name umap_self_grid_all]

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out-dir", type = "character"),
  make_option("--panel-width", type = "double", default = 4),
  make_option("--panel-height", type = "double", default = 4),
  make_option("--format", type = "character", default = "png"),
  make_option("--dpi", type = "integer", default = 300),
  make_option("--out-name", type = "character", default = "umap_self_grid_all")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

pal <- get_celltype_palette()

ORDER <- c(
  "stratamap_ill_5um_Grade1", "stratamap_ill_5um_Grade2", "stratamap_ill_5um_Grade3",
  "visiumhd", "atera", "xenium_5k", "xenium_biomarkers", "xenium_breast100"
)

# Every sample's "self" is its full QC'd gene panel ("all") - kept uniform across technologies
# (including StrataMap) for consistency with the ARI figure's SELF_GENE_SET, rather than
# restricting StrataMap to the protein-coding+lncRNA "informative" subset. Checked in the original
# pipeline that this choice barely moves downstream ARI values, and the same reasoning applies to
# the UMAP's HVG selection - the informative/all distinction mostly affects how many
# noise/technical genes (rRNA, pseudogenes, etc.) enter that computation before HVG filtering
# already narrows it down.
self_gene_set <- function(tech) "all"

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

panels <- list()
for (tech in ORDER) {
  label <- CELL_LEVEL_SAMPLES$label[CELL_LEVEL_SAMPLES$technology == tech]
  bundle_path <- reductions_path(tech, self_gene_set(tech))
  if (!file.exists(bundle_path)) { cat("SKIPPING", tech, "- no", bundle_path, "\n"); next }
  panels[[tech]] <- plot_bundle(bundle_path, label)
}

fig <- wrap_plots(panels, nrow = 1)
out_path <- file.path(opt$`out-dir`, paste0(opt$`out-name`, ".", opt$format))
width <- opt$`panel-width` * length(panels)
height <- opt$`panel-height`
if (opt$format == "png") {
  ggsave(out_path, fig, width = width, height = height, dpi = opt$dpi, limitsize = FALSE)
} else {
  ggsave(out_path, fig, width = width, height = height, limitsize = FALSE)
}
cat("Saved:", out_path, "\n")
