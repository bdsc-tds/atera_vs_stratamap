#!/usr/bin/env Rscript
# Figure 4: which cell types are most affected by raising the nCount QC floor above the pipeline's
# default (>= 10, qc_seurat.R --min-counts), for all 8 cell-level samples. Cell types with naturally
# low RNA content lose a disproportionate share of their cells as the floor goes up.
#
# --layout grid (default) produces a 2x4 facet grid; --layout row/row_tall produce 1x8 single-row
# variants at different per-panel aspect ratios. Figure 4's exact target is row_tall.
#
# Usage: Rscript celltype_threshold_sensitivity.R --extract-dir <count_distributions_raw dir> \
#   --out-dir <dir> [--layout grid|row|row_tall]

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--extract-dir", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--layout", type = "character", default = "grid")
)))
stopifnot(opt$layout %in% c("grid", "row", "row_tall"))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

THRESHOLDS <- c(10, 20, 50, 100, 200, 500, 1000, 2000)
REPORT_THRESHOLD <- 200
pal_ct <- get_celltype_palette()

all_rows <- list()
for (tech in CELL_LEVEL_SAMPLES$technology) {
  f <- file.path(opt$`extract-dir`, paste0(tech, ".parquet"))
  if (!file.exists(f)) { cat("SKIPPING", tech, "- no extract at", f, "\n"); next }
  d <- read_parquet(f) %>% filter(gene_set == "all", !is.na(first_type))
  if (nrow(d) == 0) { cat("SKIPPING", tech, "- no gene_set=='all' rows\n"); next }

  baseline_counts <- d %>% count(first_type, name = "n_total")
  retention <- lapply(THRESHOLDS, function(t) {
    d %>% filter(ncount >= t) %>% count(first_type, name = "n_retained") %>% mutate(threshold = t)
  }) %>%
    bind_rows() %>%
    right_join(baseline_counts, by = "first_type") %>%
    mutate(n_retained = tidyr::replace_na(n_retained, 0), pct_retained = 100 * n_retained / n_total,
           technology = tech, label = label_for_tech(tech))
  all_rows[[tech]] <- retention
}
combined <- do.call(rbind, all_rows)
write.csv(combined, file.path(opt$`out-dir`, "celltype_threshold_sensitivity.csv"), row.names = FALSE)

present_techs <- TECH_ORDER[TECH_ORDER %in% names(all_rows)]
sample_order <- vapply(present_techs, label_for_tech, character(1))
combined$label <- factor(combined$label, levels = sample_order)
combined$first_type <- factor(combined$first_type, levels = names(pal_ct))

font_sizes <- if (opt$layout %in% c("row", "row_tall")) {
  list(axis_text = 15, axis_title = 15, strip = 9)
} else {
  list(axis_text = 12, axis_title = 14, strip = 12)
}

p <- ggplot(combined, aes(x = threshold, y = pct_retained, color = first_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10() +
  scale_color_manual(values = pal_ct, drop = FALSE, name = "Cell type") +
  {if (opt$layout %in% c("row", "row_tall")) facet_wrap(~label, nrow = 1, scales = "free_y")
   else facet_wrap(~label, ncol = 4, scales = "free_y")} +
  labs(title = NULL, x = "nCount threshold", y = "% cells retained") +
  theme_post_bar() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = font_sizes$axis_text),
        axis.text.y = element_text(size = font_sizes$axis_text),
        axis.title = element_text(size = font_sizes$axis_title),
        strip.text = element_text(size = font_sizes$strip),
        strip.background = element_blank())
out_name <- if (opt$layout == "row_tall") "celltype_threshold_sensitivity_tall.png" else "celltype_threshold_sensitivity.png"
save_post(file.path(opt$`out-dir`, out_name), p,
          switch(opt$layout, row = "threshold_sensitivity_row", row_tall = "threshold_sensitivity_row_tall",
                 "threshold_sensitivity"))
cat("\nDone.\n")
