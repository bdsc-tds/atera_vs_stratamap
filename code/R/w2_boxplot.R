#!/usr/bin/env Rscript
# Figure 6: w2 (RCTD weight_second_type) distribution across all 8 cell-level samples, as a boxplot
# in increasing-median-score order (w2_boxplot_by_score.png). Reuses the per-cell spillover cache
# from spillover.R. Also writes a per-technology 5-number-summary CSV (raw per-cell w2 values are
# not dumped - millions of rows across 8 samples - but the boxplot itself only ever shows quartiles).
#
# Usage: Rscript w2_boxplot.R --cache-root <dir with <tech>/cell_metadata_with_spillover.parquet> --out-dir <dir>

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache-root", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

rows <- list()
for (tech in CELL_LEVEL_SAMPLES$technology) {
  f <- file.path(opt$`cache-root`, tech, "cell_metadata_with_spillover.parquet")
  if (!file.exists(f)) { cat("MISSING:", f, "- skipping\n"); next }
  d <- read_parquet(f) %>% select(any_of("weight_second_type")) %>% filter(is.finite(weight_second_type))
  d$technology <- tech
  d$label <- label_for_tech(tech)
  d$tech_family <- tech_family(tech)
  rows[[tech]] <- d
}
combined <- do.call(bind_rows, rows)
fam_pal <- TECH_FAMILY_COLORS

summary_df <- combined %>% group_by(technology, label) %>%
  summarise(n = n(), min = min(weight_second_type), p25 = quantile(weight_second_type, 0.25),
            median = median(weight_second_type), p75 = quantile(weight_second_type, 0.75),
            max = max(weight_second_type), .groups = "drop")
write.csv(summary_df, file.path(opt$`out-dir`, "w2_boxplot_summary.csv"), row.names = FALSE)

score_order <- combined %>% group_by(label) %>% summarise(m = median(weight_second_type)) %>% arrange(m) %>% pull(label)

make_boxplot <- function(df, order_levels, subtitle, base_size = 9, axis_text_size = 12, axis_title_size = 13) {
  df$label <- factor(df$label, levels = rev(order_levels))
  ggplot(df, aes(x = label, y = weight_second_type, fill = tech_family)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_fill_manual(values = fam_pal, name = "Technology") +
    labs(title = NULL, subtitle = subtitle, x = NULL, y = "secondary type weight (w2)") +
    coord_flip() +
    theme_post_bar(base_size = base_size) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = axis_text_size),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_text_size),
          axis.title.x = element_text(size = axis_title_size))
}

save_post(file.path(opt$`out-dir`, "w2_boxplot_by_score.png"),
          make_boxplot(combined, score_order, "ordered by median score",
                       base_size = 6, axis_text_size = 7, axis_title_size = 7.5), "w2_score")
cat("\nDone.\n")
