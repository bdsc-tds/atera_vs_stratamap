#!/usr/bin/env Rscript
# Figure 3: single row-grid composite of the sorted-by-abundance composition bar plots (% of that
# sample's cells, cell types reordered by OVERALL abundance summed across all 8 samples), same fixed
# presentation order as spatial_grid.R. No x-axis tick labels - bar FILL already encodes cell-type
# identity via the shared palette every figure in this bundle uses. Every panel gets the same y-axis
# upper limit and a taller per-panel aspect for legibility (figure 3's exact target).
#
# Usage: Rscript composition_grid.R --summary-csv <composition_summary.csv> --out-dir <dir>

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--summary-csv", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

pal <- get_celltype_palette()
combined <- read.csv(opt$`summary-csv`)
overall_order <- combined %>% group_by(first_type) %>% summarise(total_n = sum(n)) %>%
  arrange(desc(total_n)) %>% pull(first_type)

ORDER <- c(
  "stratamap_ill_5um_Grade1", "stratamap_ill_5um_Grade2", "stratamap_ill_5um_Grade3",
  "visiumhd", "atera", "xenium_5k", "xenium_biomarkers", "xenium_breast100"
)

build_grid <- function(y_metric, y_scale, out_name, fixed_y_max = NULL, panel_width = 2, panel_height = 3,
                        y_text_size = 12, y_title_size = 13, title_size = 13, y_title_first_only = FALSE) {
  panels <- list()
  for (tech in ORDER) {
    tab <- combined %>% filter(technology == tech)
    if (nrow(tab) == 0) { cat("SKIPPING", tech, "- not in", opt$`summary-csv`, "\n"); next }
    label <- CELL_LEVEL_SAMPLES$label[CELL_LEVEL_SAMPLES$technology == tech]
    p <- build_composition_plot(tab, label, overall_order, pal, show_x_labels = FALSE,
                                 y_metric = y_metric, y_scale = y_scale) +
      theme(plot.title = element_text(size = title_size), axis.text.y = element_text(size = y_text_size),
            axis.title.y = element_text(size = y_title_size))
    if (y_title_first_only && length(panels) > 0) p <- p + labs(y = NULL)
    if (!is.null(fixed_y_max)) p <- p + coord_cartesian(ylim = c(0, fixed_y_max))
    panels[[tech]] <- p
  }
  fig <- wrap_plots(panels, nrow = 1)
  out_path <- file.path(opt$`out-dir`, out_name)
  ggsave(out_path, fig, width = panel_width * length(panels), height = panel_height, dpi = 300, limitsize = FALSE)
  cat("Saved:", out_path, "\n")
}

build_grid("pct", "linear", "composition_grid_all_pct_linear_fixedy_tall.png", fixed_y_max = max(combined$pct),
           panel_width = 2.3, panel_height = 4.025, y_text_size = 15, y_title_size = 15,
           y_title_first_only = TRUE)
