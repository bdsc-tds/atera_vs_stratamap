#!/usr/bin/env Rscript
# Figure 7: ranked bar chart of the w2-vs-spatial-neighborhood cosine similarity, all 8 cell-level
# samples, with a genuine permutation-based significance test (not just a raw cor() call) - the null
# shuffles neighborhood_weights_second_type WITHIN the same per-sample pool the observed statistic
# is pooled over (n_perm=1000, one-sided p = mean(perm >= observed)). BH correction applied across
# the 8 samples. Also computes Pearson/Spearman (cor.test()) - CSV-only, not plotted, since figure 7
# only uses cosine; see spillover_metrics_summary.csv for the full table.
#
# Usage: Rscript spillover_metrics.R --cache-root <dir with <tech>/cell_metadata_with_spillover.parquet> \
#   --out-dir <dir> [--n-perm 1000] [--seed 42] [--reuse-summary]

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--cache-root", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--n-perm", type = "integer", default = 1000),
  make_option("--seed", type = "integer", default = 42),
  make_option("--reuse-summary", action = "store_true", default = FALSE,
              help = "skip the permutation recompute and reuse an existing spillover_metrics_summary.csv in --out-dir")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

MIN_N <- 20
TRIM_FLAG_THRESHOLD <- 0.20

cosine_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
cosine_test <- function(x, y, n_perm = 1000) {
  obs <- cosine_sim(x, y)
  perm_vals <- replicate(n_perm, cosine_sim(x, sample(y)))
  # Add-one correction (Davison & Hinkley 1997): with a finite number of permutations the true
  # p-value can never be shown to be exactly 0 - it's bounded below by 1/(n_perm+1) - so count the
  # observed statistic itself as one of the "as extreme" draws rather than reporting a naive 0/n_perm.
  list(cosine = obs, p_val = (1 + sum(perm_vals >= obs)) / (n_perm + 1))
}
trim_top1 <- function(v) pmin(v, quantile(v, 0.99, na.rm = TRUE))

summary_csv <- file.path(opt$`out-dir`, "spillover_metrics_summary.csv")
if (opt$`reuse-summary` && file.exists(summary_csv)) {
  cat("--reuse-summary: loading existing", summary_csv, "\n")
  combined <- read.csv(summary_csv, stringsAsFactors = FALSE)
  combined$label <- vapply(combined$technology, label_for_tech, character(1))
  combined$tech_family <- vapply(combined$technology, tech_family, character(1))
} else {
  rows <- list()
  for (tech in CELL_LEVEL_SAMPLES$technology) {
    cache_pq <- file.path(opt$`cache-root`, tech, "cell_metadata_with_spillover.parquet")
    if (!file.exists(cache_pq)) { cat("SKIPPING", tech, "- missing", cache_pq, "\n"); next }
    d <- read_parquet(cache_pq) %>%
      select(any_of(c("weight_second_type", "neighborhood_weights_second_type"))) %>%
      filter(is.finite(weight_second_type), is.finite(neighborhood_weights_second_type))
    if (nrow(d) < MIN_N) { cat("SKIPPING", tech, "-", nrow(d), "< MIN_N =", MIN_N, "paired observations\n"); next }
    x <- d$weight_second_type; y <- d$neighborhood_weights_second_type

    ct <- cosine_test(x, y, n_perm = opt$`n-perm`)
    cos_trim <- cosine_sim(trim_top1(x), trim_top1(y))
    pct_diff <- abs(cos_trim - ct$cosine) / abs(ct$cosine)
    if (pct_diff > TRIM_FLAG_THRESHOLD) {
      cat("FLAG:", tech, "- cosine shifts", round(100 * pct_diff, 1),
          "% after trimming top 1% (", round(ct$cosine, 4), "->", round(cos_trim, 4), ") - driven by outliers, treat with caution\n")
    }

    pear <- cor.test(x, y, method = "pearson")
    spear <- cor.test(x, y, method = "spearman", exact = FALSE)

    rows[[tech]] <- data.frame(
      technology = tech, label = label_for_tech(tech), tech_family = tech_family(tech), n = nrow(d),
      cosine = ct$cosine, cosine_p = ct$p_val,
      cosine_trimmed = cos_trim, cosine_trim_pct_diff = pct_diff,
      pearson = unname(pear$estimate), pearson_p = pear$p.value,
      spearman = unname(spear$estimate), spearman_p = spear$p.value
    )
    cat(tech, ": n=", nrow(d), " cosine=", round(ct$cosine, 3), " (p=", signif(ct$p_val, 3),
        ") pearson=", round(pear$estimate, 3), " (p=", signif(pear$p.value, 3),
        ") spearman=", round(spear$estimate, 3), " (p=", signif(spear$p.value, 3), ")\n", sep = "")
  }
  combined <- do.call(rbind, rows)
  combined$cosine_q <- p.adjust(combined$cosine_p, method = "BH")
  combined$pearson_q <- p.adjust(combined$pearson_p, method = "BH")
  combined$spearman_q <- p.adjust(combined$spearman_p, method = "BH")
  write.csv(combined, summary_csv, row.names = FALSE)
}

fam_pal <- TECH_FAMILY_COLORS

plot_metric <- function(metric, title, on_bar_labels = FALSE,
                         base_size = 9, axis_text_size = 12, axis_title_size = 13, label_size = 3) {
  d <- combined %>% arrange(.data[[metric]])
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot(d, aes(x = label, y = .data[[metric]], fill = tech_family)) +
    geom_col() +
    scale_fill_manual(values = fam_pal, name = "Technology") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = NULL, x = NULL, y = title) +
    coord_flip() +
    theme_post_bar(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = axis_text_size),
          axis.text.y = element_text(angle = 0, hjust = 1, size = axis_text_size),
          axis.title.x = element_text(size = axis_title_size))
  if (on_bar_labels) {
    text_col <- "white"
    p <- p +
      geom_text(data = d, aes(y = 0, label = label, color = I(text_col)), hjust = -0.05,
                fontface = "bold", size = label_size) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            plot.margin = margin(t = 2, r = 4, b = 2, l = 10))
  }
  p
}

save_post(file.path(opt$`out-dir`, "spillover_cosine.png"),
          plot_metric("cosine", "Cosine similarity", on_bar_labels = TRUE,
                       base_size = 6, axis_text_size = 7, axis_title_size = 7.5, label_size = 2.6),
          "spillover_cosine")
cat("\nDone.\n")
