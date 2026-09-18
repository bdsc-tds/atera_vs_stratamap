#!/usr/bin/env Rscript
# Shared styling/helpers for the figure-generation scripts in this package. source() this file,
# don't Rscript it directly. One consistent, grid-assemblable look: theme_void for spatial/UMAP,
# coord_fixed for spatial, aspect.ratio=1 for UMAP, no per-panel legends, fixed per-category sizes.
#
# Every path below is resolved from three environment variables the top-level scripts set (see
# README.md) - nothing here is hardcoded to any specific machine:
#   REPRO_ROOT    - this package's own root (contains code/ and reference/)
#   REPRO_DATA    - raw + processed per-technology data (set up by 02_preprocess.sh)
#   REPRO_RESULTS - where all CSVs/figures land

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

req_env <- function(name) {
  v <- Sys.getenv(name)
  if (v == "") stop(name, " is not set - see README.md for the 3 environment variables this package needs")
  v
}
REPRO_ROOT <- req_env("REPRO_ROOT")
REPRO_DATA <- req_env("REPRO_DATA")
REPRO_RESULTS <- req_env("REPRO_RESULTS")
REFERENCE_DIR <- file.path(REPRO_ROOT, "reference")
SPATIAL_METADATA_DIR <- file.path(REPRO_RESULTS, "spatial_metadata")

TECH_DISPLAY <- c(
  atera = "ATERA", visiumhd = "VisiumHD", xenium_biomarkers = "Xenium MM",
  xenium_breast100 = "Xenium v1", xenium_5k = "Xenium 5k"
)

label_for_tech <- function(tech, samp = NA_character_) {
  if (tech %in% names(TECH_DISPLAY)) return(TECH_DISPLAY[[tech]])
  if (grepl("^stratamap_", tech)) return(gsub("^stratamap_.*_Grade([123])$", "StrataMap G\\1", tech))
  if (!is.na(samp)) return(paste(tech, samp))
  tech
}

# The 8 cell-level samples used throughout this package (StrataMap = ill_5um segmentation only).
CELL_LEVEL_SAMPLES <- data.frame(
  technology = c("atera", "visiumhd", "xenium_biomarkers", "xenium_breast100", "xenium_5k",
                 "stratamap_ill_5um_Grade1", "stratamap_ill_5um_Grade2", "stratamap_ill_5um_Grade3"),
  is_stratamap = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
CELL_LEVEL_SAMPLES$label <- vapply(CELL_LEVEL_SAMPLES$technology, label_for_tech, character(1))

# Fixed cross-bundle technology display order, for every boxplot/ridge/grid that orders samples by
# technology rather than by a data-driven statistic.
TECH_ORDER <- c("stratamap_ill_5um_Grade1", "stratamap_ill_5um_Grade2", "stratamap_ill_5um_Grade3",
                 "visiumhd", "atera", "xenium_5k", "xenium_biomarkers", "xenium_breast100")

SAMPLE_COLORS <- c(
  "ATERA" = "#4D4D4D", "VisiumHD" = "#1B9E77",
  "Xenium MM" = "#A99BC7", "Xenium 5k" = "#8064A2", "Xenium v1" = "#D4CEE8",
  "StrataMap G1" = "#D9CBB8", "StrataMap G2" = "#A9906F", "StrataMap G3" = "#6B4F3B"
)

tech_family <- function(tech) if (grepl("^stratamap_", tech)) "StrataMap" else label_for_tech(tech)

TECH_FAMILY_COLORS <- c(
  SAMPLE_COLORS[c("ATERA", "VisiumHD", "Xenium MM", "Xenium 5k", "Xenium v1")],
  "StrataMap" = unname(SAMPLE_COLORS[["StrataMap G2"]])
)

# StrataMap's full sample directory names.
STRATAMAP_SAMPLE_FULL <- c(
  Grade1 = "DCIS_IDC_Grade1-77125049", Grade2 = "IDC_Grade2-77036971", Grade3 = "IDC_Grade3-77134058"
)

EMBEDDINGS_ROOT_MAIN <- file.path(REPRO_DATA, "processed/std_seurat_analysis")

# StrataMap sits under the same {analysis-root}/{technology}/... nesting as every 10x technology -
# "ill_5um" (segmentation resolution) stands in for "cell" (annotation level), and the grade-
# qualified sample name stands in for "breast_cancer".
reductions_path <- function(technology, gene_set) {
  if (grepl("^stratamap_", technology)) {
    grade <- regmatches(technology, regexpr("Grade[123]", technology))
    sample_full <- STRATAMAP_SAMPLE_FULL[[grade]]
    return(file.path(EMBEDDINGS_ROOT_MAIN, technology, "ill_5um/breast/breast", sample_full,
                      "embeddings", paste0("reductions_", gene_set, ".rds")))
  }
  file.path(EMBEDDINGS_ROOT_MAIN, technology, "cell/breast/breast/breast_cancer/embeddings",
            paste0("reductions_", gene_set, ".rds"))
}

qc_rds_path <- function(technology) {
  if (grepl("^stratamap_", technology)) {
    grade <- regmatches(technology, regexpr("Grade[123]", technology))
    sample_full <- STRATAMAP_SAMPLE_FULL[[grade]]
    return(file.path(EMBEDDINGS_ROOT_MAIN, technology, "ill_5um/breast/breast", sample_full,
                      paste0(sample_full, "_qc_seurat.rds")))
  }
  file.path(EMBEDDINGS_ROOT_MAIN, technology, "cell/breast/breast/breast_cancer",
            "breast_cancer_qc_seurat.rds")
}

rctd_dir_path <- function(technology) {
  root_main <- file.path(REPRO_DATA, "processed/cell_type_annotation")
  if (grepl("^stratamap_", technology)) {
    grade <- regmatches(technology, regexpr("Grade[123]", technology))
    sample_full <- STRATAMAP_SAMPLE_FULL[[grade]]
    return(file.path(root_main, technology, "ill_5um/breast/breast", sample_full,
                      "lognorm/reference_based/janesick2023/rctd_class_aware/Level2/single_cell"))
  }
  file.path(root_main, technology, "cell/breast/breast/breast_cancer",
            "lognorm/reference_based/janesick2023/rctd_class_aware/Level2/single_cell")
}

# Fixed cell-type palette, persisted so every script (and every rerun) gets identical colors - built
# once from the union of first_type across all 8 cell-level spatial_*.csv extracts. Shipped
# pre-built in reference/celltype_palette.csv for exact-color reproducibility; rebuilds itself if
# that file is missing (e.g. you added a technology whose cell types weren't in the original set).
CELLTYPE_PALETTE_PATH <- file.path(REFERENCE_DIR, "celltype_palette.csv")

build_celltype_palette <- function() {
  suppressPackageStartupMessages(library(ggsci))
  files <- list.files(SPATIAL_METADATA_DIR, pattern = "^spatial_.*_cell(_Grade[123])?\\.csv$", full.names = TRUE)
  stopifnot("No cell-level spatial_*.csv files found to build the palette from" = length(files) > 0)
  all_types <- character(0)
  for (f in files) {
    d <- read.csv(f)
    all_types <- union(all_types, unique(d$first_type))
  }
  all_types <- sort(all_types[!is.na(all_types)])
  pal <- setNames(ggsci::pal_igv()(length(all_types)), all_types)
  write.csv(data.frame(cell_type = names(pal), color = unname(pal)), CELLTYPE_PALETTE_PATH, row.names = FALSE)
  pal
}

get_celltype_palette <- function() {
  if (file.exists(CELLTYPE_PALETTE_PATH)) {
    d <- read.csv(CELLTYPE_PALETTE_PATH, colClasses = "character")
    return(setNames(d$color, d$cell_type))
  }
  build_celltype_palette()
}

# Gene panel membership, for the gene mean/var plot. Priority when a gene is in multiple panels:
# Biomarkers (Xenium MM) > Breast+100 (Xenium v1) > Xenium 5k > Atera > none.
GENE_PANEL_LEVELS <- c("Xenium MM", "Xenium v1", "Xenium 5k", "ATERA", "none")
GENE_PANEL_COLORS <- c(
  "Xenium MM" = "#E64B35", "Xenium v1" = "#4DBBD5",
  "Xenium 5k" = "#00A087", "ATERA" = "#F39B7F", "none" = "black"
)

classify_gene_panel <- function(gene_symbols) {
  biomarkers <- read.csv(file.path(REFERENCE_DIR, "biomarkers_gene_list.csv"))$gene_symbol
  breast100 <- read.csv(file.path(REFERENCE_DIR, "breast100_gene_list.csv"))$gene_symbol
  xenium5k <- read.csv(file.path(REFERENCE_DIR, "xenium5k_gene_list.csv"))$gene_symbol
  atera <- read.csv(file.path(REFERENCE_DIR, "atera_gene_list.csv"))$gene_symbol

  category <- rep("none", length(gene_symbols))
  category[gene_symbols %in% atera] <- "ATERA"
  category[gene_symbols %in% xenium5k] <- "Xenium 5k"
  category[gene_symbols %in% breast100] <- "Xenium v1"
  category[gene_symbols %in% biomarkers] <- "Xenium MM"
  factor(category, levels = GENE_PANEL_LEVELS)
}

# Scale bar for spatial plots (theme_void strips axes). Assumes the caller draws y with
# scale_y_reverse() - the bar is placed near max(y) so it lands at the visual BOTTOM after reversal.
nice_scale_length <- function(span) {
  candidates <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
  candidates[which.min(abs(candidates - span * 0.2))]
}

scale_bar_layers <- function(d, x_col = "x_centroid", y_col = "y_centroid") {
  x_range <- range(d[[x_col]], na.rm = TRUE)
  y_range <- range(d[[y_col]], na.rm = TRUE)
  bar_len <- nice_scale_length(diff(x_range))
  x0 <- x_range[1] + 0.05 * diff(x_range)
  x1 <- x0 + bar_len
  y0 <- y_range[2] - 0.05 * diff(y_range)
  label <- if (bar_len >= 1000) paste0(bar_len / 1000, " mm") else paste0(bar_len, " µm")
  list(
    annotate("segment", x = x0, xend = x1, y = y0, yend = y0, linewidth = 0.8, color = "black"),
    annotate("text", x = (x0 + x1) / 2, y = y0, label = label, vjust = 1.6, size = 2.5, color = "black")
  )
}

theme_post_spatial <- function(base_size = 9) {
  theme_void(base_size = base_size) +
    theme(legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5))
}

theme_post_umap <- function(base_size = 9) {
  theme_void(base_size = base_size) +
    theme(aspect.ratio = 1, legend.position = "none", plot.title = element_text(size = 8, hjust = 0.5))
}

theme_post_bar <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(legend.position = "none", plot.title = element_text(size = 9, hjust = 0.5),
          axis.text.x = element_text(angle = 60, hjust = 1))
}

# Fixed per-category figure size, so every figure in a category is +/- the same size for grid
# assembly. Figure 1's "bigdots"/"asp2" spatial-grid variant is passed explicitly via CLI flags in
# spatial_grid.R instead of a POST_FIGURE_SIZES entry, since it's a one-off composite, not a
# per-sample category.
POST_FIGURE_SIZES <- list(
  spatial = list(width = 4, height = 4, dpi = 300),
  umap = list(width = 4, height = 4, dpi = 300),
  composition = list(width = 5, height = 5, dpi = 300),
  distributions = list(width = 6, height = 5, dpi = 300),
  w2 = list(width = 6, height = 5, dpi = 300),
  w2_score = list(width = 8 / 2.54, height = 5 / 2.54, dpi = 300),
  spillover = list(width = 6, height = 5, dpi = 300),
  spillover_cosine = list(width = 6 / 2.54, height = 4.6 / 2.54, dpi = 300),
  threshold_sensitivity = list(width = 12, height = 9, dpi = 300),
  threshold_sensitivity_row = list(width = 2.0 * 8, height = 2.1, dpi = 300),
  threshold_sensitivity_row_tall = list(width = 2.0 * 8, height = 3.5, dpi = 300)
)

# Shared spatial-panel builder, used by both the per-sample and grid spatial scripts.
spatial_metadata_file_for_tech <- function(technology) {
  if (grepl("^stratamap_", technology)) {
    grade <- regmatches(technology, regexpr("Grade[123]", technology))
    return(file.path(SPATIAL_METADATA_DIR, paste0("spatial_", technology, "_cell_", grade, ".csv")))
  }
  file.path(SPATIAL_METADATA_DIR, paste0("spatial_", technology, "_cell.csv"))
}

# atera/xenium_biomarkers ship with x/y_centroid transposed relative to the other 6 samples.
SPATIAL_FLIP_XY <- c("atera", "xenium_biomarkers")
SPATIAL_POINT_SIZE_OVERRIDE <- c(visiumhd = 0.15)

build_spatial_plot <- function(tech, label, pal, point_size = 0.08, alpha = 0.6) {
  f <- spatial_metadata_file_for_tech(tech)
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f)
  d <- d[!is.na(d$first_type), ]
  d$first_type <- factor(d$first_type, levels = names(pal))
  if (tech %in% SPATIAL_FLIP_XY) {
    d[c("x_centroid", "y_centroid")] <- d[c("y_centroid", "x_centroid")]
  }
  ggplot(d, aes(x_centroid, y_centroid, color = first_type)) +
    geom_point(size = point_size, alpha = alpha, stroke = 0) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_y_reverse() +
    coord_fixed() +
    scale_bar_layers(d) +
    labs(title = label) +
    theme_post_spatial()
}

# Shared composition-panel builder.
build_composition_plot <- function(tab, label, level_order, pal, show_x_labels = TRUE,
                                    y_metric = "n", y_scale = "sqrt") {
  stopifnot(y_metric %in% c("n", "pct"), y_scale %in% c("sqrt", "linear"))
  tab$first_type <- factor(tab$first_type, levels = level_order)
  y_lab <- paste0(if (y_metric == "n") "Cell count" else "% of cells",
                   if (y_scale == "sqrt") " (sqrt scale)" else "")
  p <- ggplot(tab, aes(x = first_type, y = .data[[y_metric]], fill = first_type)) +
    geom_col() +
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(title = label, x = NULL, y = y_lab) +
    theme_post_bar()
  if (y_scale == "sqrt") p <- p + scale_y_sqrt()
  if (!show_x_labels) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

save_post <- function(path, plot, kind) {
  sz <- POST_FIGURE_SIZES[[kind]]
  stopifnot("Unknown figure kind - add it to POST_FIGURE_SIZES in plot_common.R" = !is.null(sz))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = sz$width, height = sz$height, dpi = sz$dpi, limitsize = FALSE)
  cat("Saved:", path, "\n")
}

# Gene-set resolution, shared by every script that needs to restrict a QCed object to a named gene
# set (embeddings, sensitivity, count-distribution extraction). Kept in one place here instead of
# duplicated per script the way the internal lab version did.
gene_list_for_set <- function(set_name, all_genes) {
  if (set_name == "all") return(all_genes)
  if (set_name %in% c("informative", "lncrna")) {
    biotypes <- arrow::read_parquet(file.path(REFERENCE_DIR, "gene_biotypes.parquet"))
    symbols <- if (set_name == "informative") biotypes$gene_symbol[biotypes$is_informative]
               else biotypes$gene_symbol[biotypes$biotype == "lncRNA"]
    return(intersect(all_genes, symbols))
  }
  path_for <- c(atera_genes = "atera_gene_list.csv", panel_genes = "breast100_gene_list.csv",
                xenium_5k_genes = "xenium5k_gene_list.csv", visiumhd_genes = "visiumhd_gene_list.csv",
                shared_panel_genes = "shared_panel_gene_list.csv")
  if (set_name %in% names(path_for)) {
    genes <- read.csv(file.path(REFERENCE_DIR, path_for[[set_name]]))$gene_symbol
    return(intersect(all_genes, genes))
  }
  stop("Unknown gene set: ", set_name)
}
