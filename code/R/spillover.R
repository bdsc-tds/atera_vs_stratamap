#!/usr/bin/env Rscript
# RNA spillover estimation: does a cell's RCTD-estimated secondary cell-type weight
# (weight_second_type) track how much of that type is physically nearby
# (neighborhood_weights_second_type, from spatial k-NN)? A strong relationship means the
# "secondary" signal looks like segmentation/diffusion spillover, not real co-expression or
# doublets. Requires the SPLIT package (https://github.com/bdsc-tds/SPLIT).
#
# Not applicable to the Chromium reference (dissociated scRNA-seq, no spatial coordinates) - this
# script is only run for technologies that have x_centroid/y_centroid.
#
# Usage: Rscript spillover.R --technology atera --qc-rds <path>_qc_seurat.rds \
#   --rctd-dir <spot_results.parquet dir> --out-dir <dir> --cache-dir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(arrow)
  library(ggplot2)
  library(optparse)
})

source(file.path(Sys.getenv("REPRO_ROOT"), "code/R/reconstruct_rctd_from_parquet.R"))

option_list <- list(
  make_option("--technology", type = "character"),
  make_option("--sample", type = "character", default = NA_character_,
              help = "optional - for multi-sample technologies (e.g. StrataMap's 3 grades), adds a 'sample' column to the output CSVs"),
  make_option("--qc-rds", type = "character"),
  make_option("--rctd-dir", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--cache-dir", type = "character"),
  make_option("--rad-pruning", type = "double", default = 15),
  make_option("--k-knn", type = "integer", default = 20)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`cache-dir`, recursive = TRUE, showWarnings = FALSE)

cat("Technology:", opt$technology, "\n")
seurat_qc <- readRDS(opt$`qc-rds`)
assay_name <- DefaultAssay(seurat_qc)

rctd <- reconstruct_rctd_from_parquet(opt$`rctd-dir`)
results_df <- rctd@results$results_df
seurat_qc <- AddMetaData(seurat_qc, results_df)
cat("spot_class table:\n"); print(table(seurat_qc$spot_class, useNA = "ifany"))

sp_coords <- as.matrix(seurat_qc@meta.data %>% select(x_centroid, y_centroid))
colnames(sp_coords) <- c("ST_1", "ST_2")
seurat_qc[["spatial"]] <- CreateDimReducObject(sp_coords, assay = assay_name, key = "ST_")

sp_nw_path <- file.path(opt$`cache-dir`, "sp_nw.rds")
if (!file.exists(sp_nw_path)) {
  sp_nw <- SPLIT::build_spatial_network(
    seurat_qc, reduction = "spatial", dims = 1:2,
    DO_prune = TRUE, rad_pruning = opt$`rad-pruning`, k_knn = opt$`k-knn`
  )
  saveRDS(sp_nw, sp_nw_path)
} else {
  sp_nw <- readRDS(sp_nw_path)
}
sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = rctd)
sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)
seurat_qc <- AddMetaData(seurat_qc, sp_neigh_df)
write_parquet(seurat_qc@meta.data, file.path(opt$`cache-dir`, "cell_metadata_with_spillover.parquet"))

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = method)
}
spillover_r2 <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  summary(lm(x ~ y))$r.squared
}

exper <- seurat_qc@meta.data %>% filter(!is.na(weight_second_type), !is.na(neighborhood_weights_second_type))

overall_stats <- exper %>% summarise(
  technology = opt$technology,
  sample = opt$sample,
  pearson = safe_cor(weight_second_type, neighborhood_weights_second_type, "pearson"),
  spearman = safe_cor(weight_second_type, neighborhood_weights_second_type, "spearman"),
  r_squared = spillover_r2(weight_second_type, neighborhood_weights_second_type),
  n_cells = n()
)
write.csv(overall_stats, file.path(opt$`out-dir`, paste0("spillover_overall_", opt$technology, ".csv")), row.names = FALSE)
cat("Overall spillover:\n"); print(overall_stats)

cor_by_first_type <- exper %>% filter(!is.na(first_type)) %>% group_by(first_type) %>% summarise(
  pearson = safe_cor(weight_second_type, neighborhood_weights_second_type, "pearson"),
  r_squared = spillover_r2(weight_second_type, neighborhood_weights_second_type),
  n_cells = n(), .groups = "drop"
) %>% arrange(desc(r_squared))
cor_by_first_type$technology <- opt$technology
cor_by_first_type$sample <- opt$sample
write.csv(cor_by_first_type, file.path(opt$`out-dir`, paste0("spillover_by_first_type_", opt$technology, ".csv")), row.names = FALSE)

cor_by_second_type <- exper %>% filter(!is.na(second_type)) %>% group_by(second_type) %>% summarise(
  pearson = safe_cor(weight_second_type, neighborhood_weights_second_type, "pearson"),
  r_squared = spillover_r2(weight_second_type, neighborhood_weights_second_type),
  n_cells = n(), .groups = "drop"
) %>% arrange(desc(r_squared))
cor_by_second_type$technology <- opt$technology
cor_by_second_type$sample <- opt$sample
write.csv(cor_by_second_type, file.path(opt$`out-dir`, paste0("spillover_by_second_type_", opt$technology, ".csv")), row.names = FALSE)

w2_by_second <- exper %>% filter(!is.na(second_type)) %>% group_by(second_type) %>% summarise(
  mean_w2 = mean(weight_second_type), median_w2 = median(weight_second_type),
  p75_w2 = quantile(weight_second_type, 0.75), n_cells = n(), .groups = "drop"
)
w2_by_second$technology <- opt$technology
w2_by_second$sample <- opt$sample
write.csv(w2_by_second, file.path(opt$`out-dir`, paste0("w2_distribution_by_second_type_", opt$technology, ".csv")), row.names = FALSE)

p_overall <- ggplot(exper, aes(x = weight_second_type, y = neighborhood_weights_second_type)) +
  geom_point(alpha = 0.2, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  labs(x = "Cell-level secondary-type weight (RCTD)", y = "Secondary-type abundance in spatial neighborhood",
       title = paste0(opt$technology, ": RNA spillover (R2=", round(overall_stats$r_squared, 3), ", n=", overall_stats$n_cells, ")")) +
  theme_classic(base_size = 14)
ggsave(file.path(opt$`out-dir`, paste0("spillover_overall_scatter_", opt$technology, ".png")), p_overall, width = 8, height = 6, dpi = 200)

cat("\nDone:", opt$`out-dir`, "\n")
