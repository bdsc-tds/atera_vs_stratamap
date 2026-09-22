#!/usr/bin/env Rscript
# Combines every technology's sensitivity_by_geneset_celltype.R output into one cross-technology
# comparison CSV (this is figure 5's --overall-csv input).
#
# Usage: Rscript sensitivity_synthesis.R --sensitivity-dir <dir with <tech>/ subdirs> --out-dir <dir>

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sensitivity-dir", type = "character"),
  make_option("--out-dir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

overall_files <- list.files(opt$`sensitivity-dir`, pattern = "^sensitivity_overall_.*\\.csv$", recursive = TRUE, full.names = TRUE)
cat("Found", length(overall_files), "overall files\n")

overall_all <- do.call(rbind, lapply(overall_files, read.csv))
write.csv(overall_all, file.path(opt$`out-dir`, "sensitivity_overall_all_technologies.csv"), row.names = FALSE)
cat("Saved combined CSV:", nrow(overall_all), "overall rows\n")
