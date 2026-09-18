#!/usr/bin/env Rscript
# StrataMap cell-level QC (parameterized script, converted from an interactive R Markdown notebook
# used internally - same logic, no plots). Builds a Seurat object from
# code/python/aggregate_stratamap_barcodes_to_cells.py's output, applies the same QC thresholds
# used for every other technology in this comparison (nCount >= 10, gene detected in >= 10 cells),
# and saves the QCed object as .rds + 10x-style .h5 + .parquet.
#
# Usage: Rscript qc_seurat_stratamap.R --sample-name DCIS_IDC_Grade1-77125049 \
#   --aggregated-dir <cell-level raw_aggregated dir> --out-dir <dir> \
#   [--min-counts 10] [--min-cells-per-gene 10]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(arrow)
  library(optparse)
  library(DropletUtils)
})

option_list <- list(
  make_option("--sample-name", type = "character"),
  make_option("--aggregated-dir", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--min-counts", type = "integer", default = 10),
  make_option("--min-cells-per-gene", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))
assay_name <- "StrataMap"

counts <- Read10X(opt$`aggregated-dir`) # gene.column = 2 (symbol) default, matches RCTD reference gene-symbol matching
# Seurat silently replaces "_" with "-" in cell names on CreateSeuratObject; sanitize here first so
# cell_meta rownames still match afterwards.
colnames(counts) <- gsub("_", "-", colnames(counts))

cell_meta <- as.data.frame(read_parquet(file.path(opt$`aggregated-dir`, "cell_metadata.parquet")))
cell_meta$barcode <- gsub("_", "-", cell_meta$barcode)
rownames(cell_meta) <- cell_meta$barcode
stopifnot(all(colnames(counts) == rownames(cell_meta)))

obj <- CreateSeuratObject(counts = counts, meta.data = cell_meta, assay = assay_name)
stopifnot("x_centroid" %in% colnames(obj@meta.data))
print(obj)

nCount_col <- paste0("nCount_", assay_name)
nFeature_col <- paste0("nFeature_", assay_name)
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")
cat(nCount_col, "summary:\n"); print(summary(obj[[nCount_col]][, 1]))

gene_cells <- rowSums(GetAssayData(obj, assay = assay_name, layer = "counts") > 0)
cat("Genes detected in 0 cells:", sum(gene_cells == 0), "\n")
cat("Genes detected in <", opt$`min-cells-per-gene`, "cells:", sum(gene_cells < opt$`min-cells-per-gene`), "\n")
cells_to_drop <- sum(obj[[nCount_col]][, 1] < opt$`min-counts`)
cat("Cells with nCount <", opt$`min-counts`, ":", cells_to_drop, "of", ncol(obj), "\n")

obj <- subset(obj, subset = !!sym(nCount_col) >= opt$`min-counts`)
genes_to_keep <- names(gene_cells[gene_cells >= opt$`min-cells-per-gene`])
obj <- obj[genes_to_keep, ]
cat("After filtering - Cells:", ncol(obj), "Genes:", nrow(obj), "\n")

sp_coords <- as.matrix(obj@meta.data %>% select(x_centroid, y_centroid))
colnames(sp_coords) <- c("ST_1", "ST_2")
obj[["spatial"]] <- CreateDimReducObject(sp_coords, assay = assay_name, key = "ST_")

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
seurat_path <- file.path(opt$`out-dir`, paste0(opt$`sample-name`, "_qc_seurat.rds"))
h5_path <- file.path(opt$`out-dir`, paste0(opt$`sample-name`, "_qc_counts.h5"))
meta_path <- file.path(opt$`out-dir`, paste0(opt$`sample-name`, "_qc_metadata.parquet"))

saveRDS(obj, file = seurat_path)
DropletUtils::write10xCounts(
  path = h5_path,
  x = GetAssayData(obj, assay = assay_name, layer = "counts"),
  barcodes = colnames(obj),
  gene.id = rownames(obj),
  type = "HDF5",
  overwrite = TRUE
)
write_parquet(obj@meta.data, meta_path)

cat("Saved:\n", seurat_path, "\n", h5_path, "\n", meta_path, "\n")
