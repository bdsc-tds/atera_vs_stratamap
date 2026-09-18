#!/usr/bin/env Rscript
# Figure 1: single row-grid composite of all 8 spatial plots, fixed presentation order (StrataMap
# Grade1-3, VisiumHD, Atera, Xenium 5k, targeted v2, targeted v1 - "v2"/"v1" = Xenium Biomarkers /
# Xenium Breast+100). Each panel keeps its own coord_fixed() aspect (real tissue proportions differ
# sample to sample) but gets equal on-canvas WIDTH via patchwork::wrap_plots' default equal-column
# layout - each panel's own scale bar already conveys its own physical scale.
#
# Usage: Rscript spatial_grid.R --out-dir <dir>
#   [--point-size 0.08] [--alpha 0.6] [--panel-width 4] [--panel-height 4]
#   [--format pdf|png] [--dpi 300] [--out-name spatial_grid_all]
#
# To reproduce the "bigdots"/taller-panel look used for a social-media crop of this figure, pass
# e.g. --format png --point-size 0.3 --alpha 0.7 --panel-width 1.5 --panel-height 3
# --out-name spatial_grid_all_bigdots_asp2 (panel-height/panel-width = 2, hence "asp2").

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out-dir", type = "character"),
  make_option("--point-size", type = "double", default = 0.08),
  make_option("--alpha", type = "double", default = 0.6),
  make_option("--panel-width", type = "double", default = 4),
  make_option("--panel-height", type = "double", default = 4),
  make_option("--format", type = "character", default = "pdf"),
  make_option("--dpi", type = "integer", default = 300),
  make_option("--out-name", type = "character", default = "spatial_grid_all")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

pal <- get_celltype_palette()

ORDER <- c(
  "stratamap_ill_5um_Grade1", "stratamap_ill_5um_Grade2", "stratamap_ill_5um_Grade3",
  "visiumhd", "atera", "xenium_5k",
  "xenium_biomarkers",  # targeted v2 (newer chemistry, multimodal staining)
  "xenium_breast100"    # targeted v1 (earliest chemistry, nucleus-expansion only)
)

panels <- list()
for (tech in ORDER) {
  label <- CELL_LEVEL_SAMPLES$label[CELL_LEVEL_SAMPLES$technology == tech]
  p <- build_spatial_plot(tech, label, pal, point_size = opt$`point-size`, alpha = opt$alpha)
  if (is.null(p)) { cat("SKIPPING", tech, "- no spatial metadata\n"); next }
  panels[[tech]] <- p
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
