#!/usr/bin/env Rscript
# N genes and N cells CSVs feeding figure 5 (summary_grid_shared_panel.R draws the actual bars
# directly from these two CSVs - no PNG rendered here): per-technology, 3 rows each -
# self / panel_genes / shared_panel_genes.
#
# n_genes comes from extract_count_distributions.R's n_genes_in_set column (constant per gene_set,
# just take the first value). n_cells comes from separation_ari_comparison.R's
# separation_ari_comparison.csv (n_cells_annotated) - already has the shared_panel_genes/
# xenium_5k/xenium_biomarkers reuse logic built in.
#
# Usage: Rscript ngenes_ncells_barplot.R --extract-dir <count_distributions_raw dir> \
#   --ari-csv <separation_ari_comparison.csv> --out-dir <dir>

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(optparse)
})
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--extract-dir", type = "character"),
  make_option("--ari-csv", type = "character"),
  make_option("--out-dir", type = "character")
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

SELF_GENE_SET <- setNames(rep("all", length(TECH_ORDER)), TECH_ORDER)
SHARED_PANEL_GENE_SET <- c(
  atera = "shared_panel_genes", visiumhd = "shared_panel_genes", xenium_breast100 = "shared_panel_genes",
  stratamap_ill_5um_Grade1 = "shared_panel_genes", stratamap_ill_5um_Grade2 = "shared_panel_genes",
  stratamap_ill_5um_Grade3 = "shared_panel_genes",
  xenium_5k = "panel_genes", xenium_biomarkers = "panel_genes"
)

ngenes_rows <- list()
for (tech in TECH_ORDER) {
  f <- file.path(opt$`extract-dir`, paste0(tech, ".parquet"))
  d <- read_parquet(f)
  get_n <- function(gs) {
    v <- d %>% filter(gene_set == gs) %>% pull(n_genes_in_set) %>% unique()
    if (length(v) == 0) NA_integer_ else v[1]
  }
  ngenes_rows[[tech]] <- data.frame(
    technology = tech, label = label_for_tech(tech),
    "Self genes" = get_n(SELF_GENE_SET[[tech]]),
    "Panel genes (filtered)" = get_n("panel_genes"),
    "Shared panel genes" = get_n(SHARED_PANEL_GENE_SET[[tech]]),
    check.names = FALSE
  )
}
ngenes <- do.call(rbind, ngenes_rows) %>%
  tidyr::pivot_longer(cols = c("Self genes", "Panel genes (filtered)", "Shared panel genes"),
                       names_to = "gene_set_type", values_to = "n_genes")

ari <- read.csv(opt$`ari-csv`, stringsAsFactors = FALSE)
ncells <- ari %>% transmute(technology, gene_set_type = ari_type, n_cells = n_cells_annotated)

write.csv(ngenes, file.path(opt$`out-dir`, "ngenes_barplot.csv"), row.names = FALSE)
write.csv(ncells, file.path(opt$`out-dir`, "ncells_barplot.csv"), row.names = FALSE)
cat("\nDone.\n")
