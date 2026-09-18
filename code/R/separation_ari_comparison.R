#!/usr/bin/env Rscript
# Two ARI-based cell-type separation metrics per technology, synthesizing across the already-computed
# separation_metric_*.csv files (separation_metric.R):
#   "self"               - each technology's own native/full gene set ("all")
#   "panel_genes"        - the shared breast+100 panel_genes list, via the *_filtered embedding
#                          (gene_set == "panel_genes_filtered")
#   "shared_panel_genes" - the single literal shared-gene-list intersection (reference/
#                          shared_panel_gene_list.csv); Xenium 5k/Xenium MM reuse their own
#                          panel_genes_filtered value instead of a dedicated computation (their own
#                          panels barely overlap the shared list to begin with).
#
# Usage: Rscript separation_ari_comparison.R --results-dir <dir with separation_metric_*.csv> --out-dir <dir>

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--results-dir", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

SELF_GENE_SET <- setNames(rep("all", length(CELL_LEVEL_SAMPLES$technology)), CELL_LEVEL_SAMPLES$technology)

SHARED_PANEL_GENE_SET <- c(
  atera = "shared_panel_genes_filtered", visiumhd = "shared_panel_genes_filtered",
  xenium_breast100 = "shared_panel_genes_filtered",
  stratamap_ill_5um_Grade1 = "shared_panel_genes_filtered",
  stratamap_ill_5um_Grade2 = "shared_panel_genes_filtered",
  stratamap_ill_5um_Grade3 = "shared_panel_genes_filtered",
  xenium_5k = "panel_genes_filtered", xenium_biomarkers = "panel_genes_filtered"
)

read_one <- function(tech, csv_name) {
  path <- file.path(opt$`results-dir`, csv_name)
  if (!file.exists(path)) { cat("MISSING:", path, "\n"); return(NULL) }
  d <- read.csv(path, stringsAsFactors = FALSE)
  # Force technology to the canonical CSV_MAP key rather than trusting the file's own "technology"
  # column - StrataMap's CSV already stores the grade-qualified name, so trusting it here would
  # double-suffix ("stratamap_ill_5um_Grade1_Grade1") and break every downstream exact match.
  d <- d %>% mutate(technology = tech)
  self_gs <- SELF_GENE_SET[[tech]]
  self_row <- d %>% filter(gene_set == self_gs) %>% mutate(ari_type = "self")
  panel_row <- d %>% filter(gene_set == "panel_genes_filtered") %>% mutate(ari_type = "panel_genes")
  shared_gs <- SHARED_PANEL_GENE_SET[[tech]]
  shared_row <- d %>% filter(gene_set == shared_gs) %>% mutate(ari_type = "shared_panel_genes")
  bind_rows(self_row, panel_row, shared_row) %>% select(technology, ari_type, n_cells_annotated, ari)
}

CSV_MAP <- c(
  atera = "separation_metric_atera_cell.csv",
  visiumhd = "separation_metric_visiumhd_cell.csv",
  xenium_5k = "separation_metric_xenium_5k_cell.csv",
  xenium_biomarkers = "separation_metric_xenium_biomarkers_cell.csv",
  xenium_breast100 = "separation_metric_xenium_breast100_cell.csv",
  stratamap_ill_5um_Grade1 = "separation_metric_stratamap_ill_5um_Grade1.csv",
  stratamap_ill_5um_Grade2 = "separation_metric_stratamap_ill_5um_Grade2.csv",
  stratamap_ill_5um_Grade3 = "separation_metric_stratamap_ill_5um_Grade3.csv"
)

combined_full <- do.call(rbind, lapply(names(CSV_MAP), function(tech) read_one(tech, CSV_MAP[[tech]])))
combined_full$label <- vapply(combined_full$technology, label_for_tech, character(1))
combined_full$label <- factor(combined_full$label, levels = vapply(TECH_ORDER, label_for_tech, character(1)))
combined_full$ari_type <- factor(combined_full$ari_type, levels = c("self", "panel_genes", "shared_panel_genes"),
                                  labels = c("Self genes", "Panel genes (filtered)", "Shared panel genes"))

write.csv(combined_full, file.path(opt$`out-dir`, "separation_ari_comparison.csv"), row.names = FALSE)
cat("\nCombined ARI table:\n")
print(combined_full %>% arrange(label, ari_type), row.names = FALSE)
cat("\nDone.\n")
