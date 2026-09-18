#!/usr/bin/env Rscript
# Per-sample cell-type composition counts feeding figure 3 (composition_grid.R draws the actual
# panels from this script's CSV output - no PNG rendered here). Reuses the spatial_*.csv extracts
# from extract_spatial_metadata.R - no new compute, no Seurat object touched.
#
# Usage: Rscript celltype_composition.R --out-dir <dir>

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(make_option("--out-dir", type = "character"))))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/plot_common.R"))

all_summaries <- list()
for (i in seq_len(nrow(CELL_LEVEL_SAMPLES))) {
  tech <- CELL_LEVEL_SAMPLES$technology[i]
  f <- spatial_metadata_file_for_tech(tech)
  if (!file.exists(f)) { cat("SKIPPING", tech, "- no spatial metadata at", f, "\n"); next }
  d <- read.csv(f)
  tab <- d %>% filter(!is.na(first_type)) %>% count(first_type, name = "n") %>%
    mutate(pct = 100 * n / sum(n), technology = tech)
  all_summaries[[tech]] <- tab
}
combined <- do.call(rbind, all_summaries)
write.csv(combined, file.path(opt$`out-dir`, "composition_summary.csv"), row.names = FALSE)

overall_order <- combined %>% group_by(first_type) %>% summarise(total_n = sum(n)) %>%
  arrange(desc(total_n)) %>% pull(first_type)
cat("Overall cell-type abundance ranking:", paste(overall_order, collapse = " > "), "\n")
cat("\nDone.\n")
