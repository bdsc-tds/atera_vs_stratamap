#!/usr/bin/env Rscript
# Cell-level QC for any Xenium-family technology (Atera, Xenium targeted biomarkers, Xenium
# breast+100, Xenium 5k) or VisiumHD: nCount >= 10, gene detected in >= 10 cells, saved as an
# .rds + .h5 pair.
#
# Two input shapes handled (see --input-type):
#   h5     - native 10x cell_feature_matrix.h5 + cells.parquet (or VisiumHD's equivalent)
#   mtxdir - 10x-style mtx/features/barcodes dir + cell_metadata.parquet (Xenium 5k's/VisiumHD's
#            converted output, see code/python/convert_*.py)
# Both ship mixed feature_types in one matrix (Gene Expression + various controls) - this script
# restricts the count matrix to real genes only.
#
# Usage:
#   Rscript qc_seurat.R --technology atera --sample breast_cancer --level cell \
#     --input-type h5 --counts-path <cell_feature_matrix.h5> --metadata-path <cells.parquet> \
#     --out-dir <out_dir> [--min-counts 10] [--min-cells-per-gene 10]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(arrow)
  library(optparse)
  library(DropletUtils)
})

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character"),
  make_option("--level", type = "character", default = "cell"),
  make_option("--input-type", type = "character"), # h5 | mtxdir
  make_option("--counts-path", type = "character"),
  make_option("--metadata-path", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--min-counts", type = "integer", default = 10),
  make_option("--min-cells-per-gene", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

assay_name <- opt$technology
cat("Technology:", opt$technology, "| Sample:", opt$sample, "| Level:", opt$level, "\n")

# --- Load counts, restrict to real gene features (mixed feature_types in one matrix otherwise) ---
select_gene_matrix <- function(counts) {
  if (is.list(counts) && !is.null(names(counts))) {
    gene_key <- names(counts)[grepl("^gene", names(counts), ignore.case = TRUE)]
    stopifnot("No 'Gene Expression'/'gene' feature type found in multi-modal 10x matrix" = length(gene_key) == 1)
    return(counts[[gene_key]])
  }
  counts
}

if (opt$`input-type` == "h5") {
  counts <- Read10X_h5(opt$`counts-path`)
} else if (opt$`input-type` == "mtxdir") {
  counts <- Read10X(opt$`counts-path`)
} else {
  stop("--input-type must be 'h5' or 'mtxdir'")
}
counts <- select_gene_matrix(counts)
# Seurat silently replaces "_" with "-" in cell names on CreateSeuratObject - sanitize here first so
# cell_meta rownames still match afterwards.
colnames(counts) <- gsub("_", "-", colnames(counts))

cell_meta <- as.data.frame(read_parquet(opt$`metadata-path`))
stopifnot("metadata must have a cell_id column" = "cell_id" %in% colnames(cell_meta))
# Some raw drops store cell_id as Arrow `binary` (raw bytes), not `string` - as.data.frame() then
# gives a list-of-raw-vectors column, and rownames<- on that silently produces garbage/duplicate
# names rather than an informative error.
if (is.list(cell_meta$cell_id)) {
  cell_meta$cell_id <- vapply(cell_meta$cell_id, function(x) rawToChar(as.raw(x)), character(1))
}
cell_meta$cell_id <- gsub("_", "-", cell_meta$cell_id)
rownames(cell_meta) <- cell_meta$cell_id
cell_meta <- cell_meta[colnames(counts), , drop = FALSE]
stopifnot(all(colnames(counts) == rownames(cell_meta)))

obj <- CreateSeuratObject(counts = counts, meta.data = cell_meta, assay = assay_name)
stopifnot(all(c("x_centroid", "y_centroid") %in% colnames(obj@meta.data)))
print(obj)

# --- QC, pre-filtering ---
nCount_col <- paste0("nCount_", assay_name)
nFeature_col <- paste0("nFeature_", assay_name)
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")
cat(nCount_col, "summary:\n"); print(summary(obj[[nCount_col]][, 1]))

gene_cells <- rowSums(GetAssayData(obj, assay = assay_name, layer = "counts") > 0)
cat("Genes detected in 0 cells:", sum(gene_cells == 0), "\n")
cat("Genes detected in <", opt$`min-cells-per-gene`, "cells:", sum(gene_cells < opt$`min-cells-per-gene`), "\n")
cells_to_drop <- sum(obj[[nCount_col]][, 1] < opt$`min-counts`)
cat("Cells with nCount <", opt$`min-counts`, ":", cells_to_drop, "of", ncol(obj), "\n")

# --- Filter ---
obj <- subset(obj, subset = !!sym(nCount_col) >= opt$`min-counts`)
genes_to_keep <- names(gene_cells[gene_cells >= opt$`min-cells-per-gene`])
obj <- obj[genes_to_keep, ]
cat("After filtering - Cells:", ncol(obj), "Genes:", nrow(obj), "\n")

# --- Spatial reduction (for visualization) ---
sp_coords <- as.matrix(obj@meta.data %>% select(x_centroid, y_centroid))
colnames(sp_coords) <- c("ST_1", "ST_2")
obj[["spatial"]] <- CreateDimReducObject(sp_coords, assay = assay_name, key = "ST_")

# --- Save (.rds + .h5) ---
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
seurat_path <- file.path(opt$`out-dir`, paste0(opt$sample, "_qc_seurat.rds"))
h5_path <- file.path(opt$`out-dir`, paste0(opt$sample, "_qc_counts.h5"))

saveRDS(obj, file = seurat_path)
DropletUtils::write10xCounts(
  path = h5_path,
  x = GetAssayData(obj, assay = assay_name, layer = "counts"),
  barcodes = colnames(obj),
  gene.id = rownames(obj),
  type = "HDF5",
  overwrite = TRUE
)

cat("Saved:\n", seurat_path, "\n", h5_path, "\n")
