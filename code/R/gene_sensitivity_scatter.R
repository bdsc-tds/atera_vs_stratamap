#!/usr/bin/env Rscript
# Figure 8: per-gene sensitivity scatterplot, one point per technology, built from
# gene_sensitivity_synthesis.R's stats table - x = shared panel genes (353) median transcripts/cell,
# y = own full panel median transcripts/cell. Also writes the exact joined per-technology data as a
# CSV (gene_sensitivity_synthesis.R's stats-csv already has every scope/metric/center combination;
# this is just the one slice actually plotted).
#
# Usage: Rscript gene_sensitivity_scatter.R [--stats-csv <path>] [--out-dir <dir>]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--stats-csv", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

stats_df <- read.csv(opt$`stats-csv`, stringsAsFactors = FALSE)
stats_df$short_label <- factor(vapply(stats_df$technology, label_for_tech, character(1)),
                                levels = names(SAMPLE_COLORS))

build_scatter <- function(d, title, x_lab, y_lab, size_lab, scale_type, same_units) {
  if (scale_type == "log10") {
    if (same_units) {
      x_min <- y_min <- min(d$x, d$y, d$x_lo, d$y_lo) * 0.5
      x_max <- y_max <- max(d$x, d$y, d$x_hi, d$y_hi) * 2
    } else {
      x_min <- min(d$x, d$x_lo) * 0.5; x_max <- max(d$x, d$x_hi) * 2
      y_min <- min(d$y, d$y_lo) * 0.5; y_max <- max(d$y, d$y_hi) * 2
    }
    x_min <- max(x_min, 1e-3)
    log_labels <- trans_format("log10", math_format(10^.x))
    x_scale <- scale_x_log10(labels = log_labels, limits = c(x_min, x_max), expand = expansion(mult = c(0.02, 0.08)))
    y_scale <- scale_y_log10(labels = log_labels, limits = c(y_min, y_max))
  } else {
    if (same_units) {
      x_max <- y_max <- max(d$x, d$y, d$x_hi, d$y_hi) * 1.05
    } else {
      x_max <- max(d$x, d$x_hi) * 1.05
      y_max <- max(d$y, d$y_hi) * 1.05
    }
    x_scale <- scale_x_continuous(labels = label_number(accuracy = NULL, big.mark = ","), limits = c(0, x_max), expand = expansion(mult = c(0.01, 0.05)))
    y_scale <- scale_y_continuous(labels = label_number(accuracy = NULL, big.mark = ","), limits = c(0, y_max))
  }
  p <- ggplot(d, aes(x = x, y = y))
  if (same_units) p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8)
  p <- p +
    geom_errorbar(aes(ymin = y_lo, ymax = y_hi, color = short_label), width = 0, linewidth = 1.8, alpha = 0.6) +
    geom_errorbarh(aes(xmin = x_lo, xmax = x_hi, color = short_label), height = 0, linewidth = 1.8, alpha = 0.6) +
    geom_point(aes(size = size, fill = short_label), shape = 21, color = "black", stroke = 1, alpha = 0.85) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_size_continuous(name = size_lab, trans = "log10", labels = comma, range = c(5, 20)) +
    x_scale + y_scale +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_bw(base_size = 34) +
    theme(legend.position = "bottom", aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, face = "bold", size = 34, lineheight = 1.2),
          axis.title = element_text(size = 36, lineheight = 1.15), axis.text = element_text(size = 38),
          legend.title = element_text(size = 26), legend.text = element_text(size = 24),
          plot.margin = margin(20, 40, 20, 20))
  p
}

save_scatter <- function(p, name) {
  out_png <- file.path(opt$`out-dir`, paste0(name, ".png"))
  ggsave(out_png, p, width = 15, height = 16, dpi = 300, limitsize = FALSE)
  cat("Saved:", out_png, "\n")
}

# x = shared_panel_genes (353) median transcripts/cell, y = own_panel (self) median transcripts/cell.
own <- stats_df %>% filter(scope == "own_panel", metric == "transcripts_per_cell") %>%
  select(technology, short_label, y = median, y_lo = p25, y_hi = p75, size = n)
shared <- stats_df %>% filter(scope == "shared_panel_genes", metric == "transcripts_per_cell") %>%
  select(technology, x = median, x_lo = p25, x_hi = p75)
d <- inner_join(own, shared, by = "technology") %>% arrange(desc(size))
write.csv(d, file.path(opt$`out-dir`, "scatter_self_vs_shared_panel_genes_transcripts_per_cell_median.csv"), row.names = FALSE)

title <- "Transcripts / cell / gene (median):\nself-panel vs. shared panel genes (353)"
p <- build_scatter(d, title, "Shared panel genes (353)", "Self-panel", "Genes in self panel", "log10", same_units = TRUE)
save_scatter(p, "scatter_self_vs_shared_panel_genes_transcripts_per_cell_median_log10")
