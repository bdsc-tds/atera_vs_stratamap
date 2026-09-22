#!/usr/bin/env Rscript
# Synthesizes gene_sensitivity.R's per-gene, per-technology output into summary statistics
# (n/mean/median/sd/p25/p75/min/max) for two gene scopes:
#   - "own_panel": every gene in the technology's own post-QC panel (raw output, unchanged).
#   - "shared_panel_genes": restricted to the literal 353-gene shared-gene-list intersection, the
#     SAME list and the SAME filter for every technology - including Xenium 5k/Xenium MM, which
#     previously reused a proxy (their own panel intersected with the 380-gene breast+100 list
#     instead) that both understated the true overlap and, worse, wasn't the quantity figure 8's
#     axis label ("Shared panel genes (353)") claimed it was. Xenium 5k actually has 229 of the 353
#     genes in its own panel, Xenium MM has 91 - both plenty for a real median/IQR, so there was no
#     need for a proxy at all.
#
# Usage: Rscript gene_sensitivity_synthesis.R --sensitivity-dir <gene_sensitivity.R's --out-dir> --out-dir <dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sensitivity-dir", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(file.path(opt$`out-dir`, "csv"), recursive = TRUE, showWarnings = FALSE)

combined <- read.csv(file.path(opt$`sensitivity-dir`, "gene_sensitivity_combined.csv"))
shared_panel_genes <- read.csv(file.path(REFERENCE_DIR, "shared_panel_gene_list.csv"))$gene_symbol
cat("Loaded", nrow(combined), "gene x technology rows across", length(unique(combined$technology)), "technologies\n")

scopes <- list(
  own_panel = combined,
  shared_panel_genes = combined %>% filter(gene %in% shared_panel_genes)
)
for (s in names(scopes)) cat("Scope", s, ":", nrow(scopes[[s]]), "rows,",
                              length(unique(scopes[[s]]$technology)), "technologies with >0 rows\n")

metrics <- c("transcripts_per_cell", "transcripts_per_mm2")

stats_all <- list()
for (scope_name in names(scopes)) {
  df <- scopes[[scope_name]]
  for (metric in metrics) {
    st <- df %>%
      group_by(technology, label) %>%
      summarise(
        n = n(), mean = mean(.data[[metric]]), median = median(.data[[metric]]), sd = sd(.data[[metric]]),
        p25 = quantile(.data[[metric]], 0.25), p75 = quantile(.data[[metric]], 0.75),
        min = min(.data[[metric]]), max = max(.data[[metric]]), .groups = "drop"
      ) %>%
      mutate(scope = scope_name, metric = metric)
    stats_all[[paste(scope_name, metric)]] <- st
  }
}
stats_df <- bind_rows(stats_all) %>%
  select(scope, metric, technology, label, n, mean, median, sd, p25, p75, min, max)
out_csv <- file.path(opt$`out-dir`, "csv", "gene_sensitivity_stats.csv")
write.csv(stats_df, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")
