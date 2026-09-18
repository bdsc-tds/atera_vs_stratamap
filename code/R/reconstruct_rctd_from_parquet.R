# Reconstructs a minimal `RCTD` S4 object (from the `spacexr` package) out of the parquet/h5 files
# written by code/python/run_rctd_doublet.py, so downstream R code (spillover.R) can use SPLIT's
# functions, which expect a real RCTD object rather than plain data frames. Only the fields SPLIT
# actually reads are populated (spatial coordinates are dummied out at 0,0 - not used by
# SPLIT::build_spatial_network here, which is called with the real coordinates from a separate
# Seurat dimensionality reduction instead).
#
# Requires: spacexr (https://github.com/dmcable/spacexr) and SPLIT (https://github.com/bdsc-tds/SPLIT)
reconstruct_rctd_from_parquet <- function(save_dir) {
  library(arrow)
  library(rhdf5)
  library(Matrix)
  library(spacexr)
  library(SPLIT)

  cell_ids_df     <- read_parquet(file.path(save_dir, "cell_ids.parquet"))
  weights_df      <- read_parquet(file.path(save_dir, "weights.parquet"))
  weights_dbl_df  <- read_parquet(file.path(save_dir, "weights_doublet.parquet"))
  spot_results_df <- read_parquet(file.path(save_dir, "spot_results.parquet"))
  pixel_mask_df   <- read_parquet(file.path(save_dir, "pixel_mask.parquet"))
  metadata_df     <- read_parquet(file.path(save_dir, "metadata.parquet"))

  cell_ids        <- cell_ids_df$cell_id
  cell_type_names <- metadata_df$cell_type_names

  spot_class_map <- c("0" = "reject", "1" = "singlet",
                      "2" = "doublet_certain", "3" = "doublet_uncertain")

  results_df <- data.frame(
    spot_class    = factor(
      spot_class_map[as.character(spot_results_df$spot_class)],
      levels = c("singlet", "doublet_certain", "doublet_uncertain", "reject")
    ),
    first_type    = factor(spot_results_df$first_type_name,  levels = cell_type_names),
    second_type   = factor(spot_results_df$second_type_name, levels = cell_type_names),
    first_class   = factor(
      ifelse(spot_results_df$first_class, "doublet_certain", "singlet"),
      levels = c("singlet", "doublet_certain", "doublet_uncertain", "reject")
    ),
    second_class  = factor(
      ifelse(spot_results_df$second_class, "doublet_certain", "singlet"),
      levels = c("singlet", "doublet_certain", "doublet_uncertain", "reject")
    ),
    min_score     = spot_results_df$min_score,
    singlet_score = spot_results_df$singlet_score,
    row.names     = cell_ids
  )
  stopifnot(!any(is.na(results_df$spot_class)))

  weights_mat <- as.matrix(weights_df)
  rownames(weights_mat) <- cell_ids
  colnames(weights_mat) <- cell_type_names

  weights_doublet_mat <- as.matrix(weights_dbl_df[, c("w_1", "w_2")])
  rownames(weights_doublet_mat) <- cell_ids
  colnames(weights_doublet_mat) <- c("first_type", "second_type")

  h5_path        <- file.path(save_dir, "reference_profiles.h5")
  prof           <- h5read(h5_path, "profiles")
  gene_names_ref <- h5read(h5_path, "gene_names")
  ct_names_ref   <- h5read(h5_path, "cell_type_names")

  # prof is cell_types x genes from Python, transpose to genes x cell_types
  profiles_df <- as.data.frame(t(prof))
  rownames(profiles_df) <- gene_names_ref
  colnames(profiles_df) <- ct_names_ref

  cell_type_info <- list(
    info   = list(profiles_df, ct_names_ref, length(ct_names_ref)),
    renorm = NULL
  )

  all_cell_ids <- pixel_mask_df$cell_id
  n_cells      <- length(all_cell_ids)

  coords <- data.frame(
    x = rep(0, n_cells), y = rep(0, n_cells),
    row.names = all_cell_ids
  )
  dummy_counts <- Matrix::Matrix(
    0L, nrow = 1, ncol = n_cells,
    dimnames = list("dummy_gene", all_cell_ids),
    sparse = TRUE
  )
  spatialRNA_obj <- new(
    "SpatialRNA",
    coords = coords,
    counts = dummy_counts,
    nUMI   = setNames(rep(1L, n_cells), all_cell_ids)
  )

  RCTD <- new(
    "RCTD",
    spatialRNA         = spatialRNA_obj,
    originalSpatialRNA = spatialRNA_obj,
    reference          = new("Reference"),
    results            = list(
      results_df_xe   = results_df,
      weights         = weights_mat,
      weights_doublet = weights_doublet_mat
    ),
    cell_type_info = cell_type_info,
    config         = list(RCTDmode = "doublet"),
    internal_vars  = list(
      class_df = data.frame(class = cell_type_names, row.names = cell_type_names)
    ),
    de_results     = list(),
    internal_vars_de = list()
  )

  SPLIT::run_post_process_RCTD(RCTD, lite = TRUE, min_weight = 0.01)
}
