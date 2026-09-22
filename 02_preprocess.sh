#!/bin/bash
# Raw data -> every processed artifact the 8 figures in 03_generate_figures.sh depend on: per-
# technology ingestion -> QC -> RCTD cell-type annotation -> gene-set embeddings -> spatial
# metadata -> spillover cache -> composition/sensitivity/ARI summaries, for all 6 technologies
# (Atera, VisiumHD, Xenium MM/biomarkers, Xenium v1/breast+100, Xenium 5k, StrataMap x3 grades).
#
# Run as: bash 02_preprocess.sh (not sourced). Safe to re-run - every step is idempotent (same
# output paths) and most scripts just overwrite their outputs.
#
# Prerequisites: an R + Python environment with the packages in environment.yml (conda env create
# -f environment.yml && conda activate spatial-tech-comparison), plus 2 GitHub-only R packages not
# on CRAN/Bioconductor (install once):
#   R -e 'remotes::install_github("dmcable/spacexr")'
#   R -e 'remotes::install_github("bdsc-tds/SPLIT")'
# See 00_download_raw_data.md for where every raw dataset comes from.

set -euo pipefail

# ============================================================ CONFIG
# REPRO_ROOT/REPRO_DATA/REPRO_RESULTS, USE_SLURM/SLURM_ACCOUNT/SLURM_PARTITION, R_CMD/PYTHON_CMD,
# and the dispatch()/wait_all() helpers all live in code/lib.sh (shared with
# 03_generate_figures.sh) - see that file for what each one does and how to override it.
# REPRO_DATA will be LARGE (raw downloads alone run tens of GB; StrataMap's raw barcode matrices
# are the biggest single piece) - point it at a disk with room, not necessarily inside REPRO_ROOT.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/code/lib.sh"

RAW="$REPRO_DATA/raw"
PROC="$REPRO_DATA/processed"
STD_ROOT="$PROC/std_seurat_analysis"
RCTD_ROOT="$PROC/cell_type_annotation"
SEG_ROOT="$PROC/segmentation"
REFERENCE_H5AD="$REPRO_DATA/reference/chromium_flex_reference.h5ad"
mkdir -p "$RAW" "$PROC"

declare -A SM_SAMPLE_FULL=( [Grade1]=DCIS_IDC_Grade1-77125049 [Grade2]=IDC_Grade2-77036971 [Grade3]=IDC_Grade3-77134058 )

# StrataMap uses the same {analysis-root}/{technology}/... nesting as every 10x technology below -
# just with "ill_5um" (segmentation resolution) standing in for "cell" (annotation level) and the
# grade-qualified sample name standing in for "breast_cancer".
qc_rds_for()  { local t=$1; [[ $t == stratamap_* ]] && { local g=${t##*_}; echo "$STD_ROOT/$t/ill_5um/breast/breast/${SM_SAMPLE_FULL[$g]}/${SM_SAMPLE_FULL[$g]}_qc_seurat.rds"; return; }; echo "$STD_ROOT/$t/cell/breast/breast/breast_cancer/breast_cancer_qc_seurat.rds"; }
rctd_dir_for(){ local t=$1; [[ $t == stratamap_* ]] && { local g=${t##*_}; echo "$RCTD_ROOT/$t/ill_5um/breast/breast/${SM_SAMPLE_FULL[$g]}/lognorm/reference_based/janesick2023/rctd_class_aware/Level2/single_cell"; return; }; echo "$RCTD_ROOT/$t/cell/breast/breast/breast_cancer/lognorm/reference_based/janesick2023/rctd_class_aware/Level2/single_cell"; }
embed_dir_for(){ local t=$1; [[ $t == stratamap_* ]] && { local g=${t##*_}; echo "$STD_ROOT/$t/ill_5um/breast/breast/${SM_SAMPLE_FULL[$g]}/embeddings"; return; }; echo "$STD_ROOT/$t/cell/breast/breast/breast_cancer/embeddings"; }
seg_dir_for() { local t=$1; local g=${t##*_}; echo "$SEG_ROOT/$t/ill_5um/breast/breast/${SM_SAMPLE_FULL[$g]}/raw_aggregated"; }

TECHS=(atera xenium_biomarkers xenium_breast100 xenium_5k visiumhd)
SAMPLES=(atera xenium_biomarkers xenium_breast100 xenium_5k visiumhd stratamap_ill_5um_Grade1 stratamap_ill_5um_Grade2 stratamap_ill_5um_Grade3)

# ============================================================ STAGE 0: Chromium reference +
# VisiumHD download (the two pieces this script can fetch itself with no bot-check/manual step).
echo "=== STAGE 0: Chromium reference + VisiumHD download ==="
if [ ! -f "$REFERENCE_H5AD" ]; then
  dispatch build_reference 16G 00:30:00 2 \
    "$PYTHON_CMD $CODE_PY/build_reference.py --out $REFERENCE_H5AD --geo-h5 $RAW/reference/GSM7782698_count_raw_feature_bc_matrix.h5 --celltype-xlsx $RAW/reference/Cell_Barcode_Type_Matrices.xlsx --taxonomy-csv $REPRO_ROOT/reference/janesick_celltype_taxonomy.csv"
  wait_all
fi

mkdir -p "$RAW/visiumhd"
VHD_TGZ="$RAW/visiumhd/Visium_HD_11mm_Human_Breast_Cancer_segmented_outputs.tar.gz"
if [ ! -f "$VHD_TGZ" ]; then
  curl -o "$VHD_TGZ" "https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_11mm_Human_Breast_Cancer/Visium_HD_11mm_Human_Breast_Cancer_segmented_outputs.tar.gz"
fi
mkdir -p "$RAW/visiumhd/segmented_outputs_extracted"
VISIUMHD_SEG_EXTRACTED="$RAW/visiumhd/segmented_outputs_extracted/segmented_outputs"
if [ ! -f "$VISIUMHD_SEG_EXTRACTED/filtered_feature_cell_matrix.h5" ]; then
  tar -xzf "$VHD_TGZ" -C "$RAW/visiumhd/segmented_outputs_extracted"
fi

# The other 4 raw drops (Atera, Xenium MM, Xenium v1, Xenium 5k) and StrataMap need a manual
# browser download - see 00_download_raw_data.md. This script assumes they're already unpacked
# under $REPRO_DATA/raw/ by the time you reach STAGE 1.

# ============================================================ STAGE 1: per-technology ingestion
# (Xenium 5k's zarr bundle, VisiumHD's 2um-bin-to-cell aggregation, StrataMap's raw-barcode-to-cell
# aggregation). Feeds STAGE 2.
echo "=== STAGE 1: ingestion ==="
XENIUM5K_CONVERTED="$SEG_ROOT/xenium_5k/cell/breast/breast/breast_cancer"
VISIUMHD_CONVERTED="$SEG_ROOT/visiumhd/cell/breast/breast/breast_cancer"

dispatch xe5k_convert 32G 00:30:00 4 \
  "$PYTHON_CMD $CODE_PY/convert_xenium5k_zarr_cells.py $RAW/xenium_5k/xe_outs $XENIUM5K_CONVERTED"
# VisiumHD's geojson centroids are in full-res image pixels, not microns - the converter rescales
# using spatial/scalefactors_json.json's microns_per_pixel (required 3rd positional arg).
dispatch vhd_convert 32G 01:00:00 4 \
  "$PYTHON_CMD $CODE_PY/convert_visiumhd_segmented_cells.py $VISIUMHD_SEG_EXTRACTED $VISIUMHD_CONVERTED $VISIUMHD_SEG_EXTRACTED/spatial/scalefactors_json.json"
wait_all

declare -A SM_AGG_MEM=( [Grade1]=32G [Grade2]=64G [Grade3]=250G )
for GRADE in Grade1 Grade2 Grade3; do
  SAMPLE_FULL=${SM_SAMPLE_FULL[$GRADE]}
  dispatch "sm_aggregate_$GRADE" "${SM_AGG_MEM[$GRADE]}" 04:00:00 4 \
    "$PYTHON_CMD $CODE_PY/aggregate_stratamap_barcodes_to_cells.py \
      --raw-dir $RAW/stratamap/$SAMPLE_FULL/Raw_Matrix_Files \
      --segmentation-csv $RAW/stratamap/$SAMPLE_FULL/Segmentations/registered_Expanded_5um_cell_contour_coords_local.csv \
      --out-dir $(seg_dir_for stratamap_ill_5um_$GRADE)"
done
wait_all

declare -A COUNTS_PATH=(
  [atera]="$RAW/atera/outs/cell_feature_matrix.h5"
  [xenium_biomarkers]="$RAW/xenium_biomarkers/outs/cell_feature_matrix.h5"
  [xenium_breast100]="$RAW/xenium_breast100/outs/cell_feature_matrix.h5"
  [xenium_5k]="$XENIUM5K_CONVERTED"
  [visiumhd]="$VISIUMHD_SEG_EXTRACTED/filtered_feature_cell_matrix.h5"
)
declare -A METADATA_PATH=(
  [atera]="$RAW/atera/outs/cells.parquet"
  [xenium_biomarkers]="$RAW/xenium_biomarkers/outs/cells.parquet"
  [xenium_breast100]="$RAW/xenium_breast100/outs/cells.parquet"
  [xenium_5k]="$XENIUM5K_CONVERTED/cell_metadata.parquet"
  [visiumhd]="$VISIUMHD_CONVERTED/cell_metadata.parquet"
)
declare -A INPUT_TYPE=( [atera]=h5 [xenium_biomarkers]=h5 [xenium_breast100]=h5 [xenium_5k]=mtxdir [visiumhd]=h5 )

# ============================================================ STAGE 2: QC
echo "=== STAGE 2: QC ==="
declare -A QC_MEM=( [atera]=32G [xenium_biomarkers]=32G [xenium_breast100]=32G [xenium_5k]=48G [visiumhd]=64G )
for TECH in "${TECHS[@]}"; do
  OUT_DIR="$STD_ROOT/$TECH/cell/breast/breast/breast_cancer"
  dispatch "qc_$TECH" "${QC_MEM[$TECH]}" 01:00:00 4 \
    "$R_CMD $CODE_R/qc_seurat.R --technology $TECH --sample breast_cancer --level cell \
      --input-type ${INPUT_TYPE[$TECH]} --counts-path ${COUNTS_PATH[$TECH]} --metadata-path ${METADATA_PATH[$TECH]} \
      --out-dir $OUT_DIR"
done
wait_all

declare -A SM_QC_MEM=( [Grade1]=64G [Grade2]=96G [Grade3]=150G )
for GRADE in Grade1 Grade2 Grade3; do
  SAMPLE_FULL=${SM_SAMPLE_FULL[$GRADE]}
  dispatch "qc_sm_$GRADE" "${SM_QC_MEM[$GRADE]}" 01:30:00 4 \
    "$R_CMD $CODE_R/qc_seurat_stratamap.R --sample-name $SAMPLE_FULL \
      --aggregated-dir $(seg_dir_for stratamap_ill_5um_$GRADE) \
      --out-dir $(dirname "$(qc_rds_for stratamap_ill_5um_$GRADE)")"
done
wait_all

# ============================================================ STAGE 3: RCTD
echo "=== STAGE 3: RCTD ==="
declare -A RCTD_MEM=( [atera]=64G [xenium_biomarkers]=64G [xenium_breast100]=64G [xenium_5k]=64G [visiumhd]=250G )
for TECH in "${TECHS[@]}"; do
  QC_COUNTS_H5="$STD_ROOT/$TECH/cell/breast/breast/breast_cancer/breast_cancer_qc_counts.h5"
  dispatch "rctd_$TECH" "${RCTD_MEM[$TECH]}" 02:00:00 8 \
    "$PYTHON_CMD $CODE_PY/run_rctd_doublet.py --h5 $QC_COUNTS_H5 --reference-h5ad $REFERENCE_H5AD --out-dir $(rctd_dir_for $TECH)"
done
wait_all

declare -A SM_RCTD_MEM=( [Grade1]=64G [Grade2]=96G [Grade3]=250G )
for GRADE in Grade1 Grade2 Grade3; do
  SAMPLE_FULL=${SM_SAMPLE_FULL[$GRADE]}
  QC_COUNTS_H5="$(dirname "$(qc_rds_for stratamap_ill_5um_$GRADE)")/${SAMPLE_FULL}_qc_counts.h5"
  dispatch "rctd_sm_$GRADE" "${SM_RCTD_MEM[$GRADE]}" 04:00:00 8 \
    "$PYTHON_CMD $CODE_PY/run_rctd_doublet.py --h5 $QC_COUNTS_H5 --reference-h5ad $REFERENCE_H5AD --out-dir $(rctd_dir_for stratamap_ill_5um_$GRADE)"
done
wait_all

# ============================================================ STAGE 4: gene-set embeddings.
# Feeds fig 2 (umap_self_*) and, via separation metric/ARI below, fig 5. Trimmed to just the "self"
# gene set - atera_genes/panel_genes embeddings were computed by the original internal pipeline but
# neither figure 2 (only uses the "self" UMAP group) nor figure 5 (only ever plots "Self genes" vs.
# "Shared panel genes" ARI/n_genes/n_cells bars - checked summary_grid_shared_panel.R's draw_panel_*
# functions directly, panel_genes/panel_genes_filtered are never read) actually consumes them - real
# dead compute for this figure set, not just "nice to trim".
echo "=== STAGE 4: gene-set embeddings (self only) ==="
declare -A EMBED_MEM=( [atera]=64G [xenium_biomarkers]=64G [xenium_breast100]=64G [xenium_5k]=64G [visiumhd]=96G )
for TECH in "${TECHS[@]}"; do
  # (umap_self_grid.R's self_gene_set()) - Atera/VisiumHD's self group uses "all" too, so no
  # technology in this loop needs it.
  dispatch "embed_$TECH" "${EMBED_MEM[$TECH]}" 01:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample breast_cancer \
      --qc-rds $(qc_rds_for $TECH) --gene-sets all --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH)"
done
wait_all

declare -A SM_EMBED_MEM=( [Grade1]=64G [Grade2]=110G [Grade3]=220G )
for GRADE in Grade1 Grade2 Grade3; do
  TECH="stratamap_ill_5um_$GRADE"
  UMAP_INIT=""; [ "$GRADE" = "Grade3" ] && UMAP_INIT="--umap-init random"  # RSpectra::eigs_sym can segfault on ~700k-cell inputs with the spectral default
  # StrataMap's "self" UMAP uses "all", matching every other
  # technology and the ARI figure's SELF_GENE_SET - see umap_self_grid.R's self_gene_set().
  dispatch "embed_sm_$GRADE" "${SM_EMBED_MEM[$GRADE]}" 02:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample $GRADE \
      --qc-rds $(qc_rds_for $TECH) --gene-sets all \
      --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH) $UMAP_INIT"
done
wait_all

# panel_genes_filtered - ONLY needed for xenium_5k/xenium_biomarkers, as their proxy for the
# "shared_panel_genes" ARI bar (they don't get a native shared_panel_genes computation below since
# their own panels barely overlap the breast+100 gene list to begin with - see STAGE 4c). Every
# other technology gets shared_panel_genes_filtered directly in STAGE 4c instead, so this step
# shrinks from "all 8 samples" to just these 2.
echo "=== STAGE 4b: panel_genes_filtered embeddings (xenium_5k, xenium_biomarkers only - shared-panel proxy) ==="
for TECH in xenium_5k xenium_biomarkers; do
  # Xenium MM's own panel is only 96 genes, below gene_set_embeddings.R's --min-genes default of
  # 200 - it would otherwise be silently skipped even though 96 genes is plenty for a PCA/ARI here.
  MIN_GENES=200; [ "$TECH" = "xenium_biomarkers" ] && MIN_GENES=50
  dispatch "embed_filt_$TECH" 64G 01:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample breast_cancer \
      --qc-rds $(qc_rds_for $TECH) --gene-sets panel_genes --min-count-in-subset 10 --npcs 15 --label-suffix _filtered \
      --min-genes $MIN_GENES --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH)"
done
wait_all

# shared_panel_genes_filtered - only for the technologies with a real (non-proxy) shared axis;
# xenium_5k/xenium_biomarkers reuse their own panel_genes_filtered as this proxy instead.
# --min-count-in-subset 10 matches the >=10 nCount QC floor already applied to the self/"all" panel
# (qc_seurat.R's --min-counts default) - the same cell-inclusion threshold everywhere, not a
# different, stricter one just for the shared-panel-genes axis.
echo "=== STAGE 4c: shared_panel_genes_filtered embeddings (6 samples) ==="
for TECH in atera xenium_breast100 visiumhd stratamap_ill_5um_Grade1 stratamap_ill_5um_Grade2 stratamap_ill_5um_Grade3; do
  MEM=64G; [[ $TECH == stratamap_*Grade2 ]] && MEM=110G; [[ $TECH == stratamap_*Grade3 ]] && MEM=220G
  UMAP_INIT=""; [[ $TECH == stratamap_*Grade3 ]] && UMAP_INIT="--umap-init random"
  SAMPLE=breast_cancer; [[ $TECH == stratamap_* ]] && SAMPLE=${TECH##*_}
  dispatch "embed_shared_$TECH" "$MEM" 02:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample $SAMPLE \
      --qc-rds $(qc_rds_for $TECH) --gene-sets shared_panel_genes --min-count-in-subset 10 --npcs 15 --label-suffix _filtered \
      --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH) $UMAP_INIT"
done
wait_all

# ============================================================ STAGE 5: spatial metadata extract.
# Feeds fig 1 directly, fig 8 indirectly (gene_sensitivity.R requires this to exist first).
echo "=== STAGE 5: spatial metadata ==="
for TECH in "${TECHS[@]}"; do
  dispatch "spatial_$TECH" 48G 00:45:00 2 \
    "$R_CMD $CODE_R/extract_spatial_metadata.R --technology $TECH --level cell \
      --qc-rds $(qc_rds_for $TECH) --rctd-dir $(rctd_dir_for $TECH) --out-dir $REPRO_RESULTS/spatial_metadata"
done
for GRADE in Grade1 Grade2 Grade3; do
  TECH="stratamap_ill_5um_$GRADE"
  dispatch "spatial_sm_$GRADE" 64G 00:45:00 2 \
    "$R_CMD $CODE_R/extract_spatial_metadata.R --technology $TECH --sample $GRADE --level cell \
      --qc-rds $(qc_rds_for $TECH) --rctd-dir $(rctd_dir_for $TECH) --out-dir $REPRO_RESULTS/spatial_metadata"
done
wait_all

# ============================================================ STAGE 6: spillover cache. Feeds
# figs 6/7 directly.
echo "=== STAGE 6: spillover ==="
declare -A SPILL_MEM=( [atera]=64G [xenium_biomarkers]=64G [xenium_breast100]=64G [xenium_5k]=64G [visiumhd]=96G )
for TECH in "${TECHS[@]}"; do
  dispatch "spillover_$TECH" "${SPILL_MEM[$TECH]}" 03:00:00 8 \
    "$R_CMD $CODE_R/spillover.R --technology $TECH \
      --qc-rds $(qc_rds_for $TECH) --rctd-dir $(rctd_dir_for $TECH) \
      --cache-dir $REPRO_RESULTS/spillover_cache/$TECH"
done
declare -A SM_SPILL_MEM=( [Grade1]=64G [Grade2]=96G [Grade3]=150G )
for GRADE in Grade1 Grade2 Grade3; do
  # Cache dir naming below (stratamap_ill_5um_<Grade>) matches plot_common.R's tech key exactly.
  dispatch "spillover_sm_$GRADE" "${SM_SPILL_MEM[$GRADE]}" 04:00:00 8 \
    "$R_CMD $CODE_R/spillover.R --technology stratamap_ill_5um_$GRADE \
      --qc-rds $(qc_rds_for stratamap_ill_5um_$GRADE) --rctd-dir $(rctd_dir_for stratamap_ill_5um_$GRADE) \
      --cache-dir $REPRO_RESULTS/spillover_cache/stratamap_ill_5um_$GRADE"
done
wait_all

# ============================================================ STAGE 7: composition summary.
# Feeds fig 3.
echo "=== STAGE 7: composition summary ==="
dispatch composition 16G 00:15:00 2 \
  "$R_CMD $CODE_R/celltype_composition.R --out-dir $REPRO_RESULTS/composition"
wait_all

# ============================================================ STAGE 8: separation metric + ARI
# comparison. Feeds fig 5's --ari-csv. separation_metric.R does real compute (k-means on up to
# ~700k-cell PCA embeddings) - not a "cheap plotting" step, dispatched like everything else.
echo "=== STAGE 8: separation metric (ARI) ==="
for TECH in "${TECHS[@]}"; do
  dispatch "sepmetric_$TECH" 32G 00:30:00 4 \
    "$R_CMD $CODE_R/separation_metric.R --technology $TECH \
      --embeddings-dir $(embed_dir_for "$TECH") --out-csv $REPRO_RESULTS/separation_metric/separation_metric_${TECH}_cell.csv"
done
for GRADE in Grade1 Grade2 Grade3; do
  dispatch "sepmetric_sm_$GRADE" 32G 00:30:00 4 \
    "$R_CMD $CODE_R/separation_metric.R --technology stratamap_ill_5um_$GRADE \
      --embeddings-dir $(embed_dir_for stratamap_ill_5um_$GRADE) --out-csv $REPRO_RESULTS/separation_metric/separation_metric_stratamap_ill_5um_$GRADE.csv"
done
wait_all
dispatch separation_ari 16G 00:15:00 2 \
  "$R_CMD $CODE_R/separation_ari_comparison.R --results-dir $REPRO_RESULTS/separation_metric --out-dir $REPRO_RESULTS/separation_ari"
wait_all

# ============================================================ STAGE 9: per-cell count-distribution
# extraction. Feeds fig 4 directly, fig 5's ngenes/ncells panels (STAGE 10).
echo "=== STAGE 9: count distribution extraction ==="
declare -A EXTRACT_MEM=( [atera]=48G [xenium_biomarkers]=48G [xenium_breast100]=48G [xenium_5k]=48G [visiumhd]=64G \
  [stratamap_ill_5um_Grade1]=64G [stratamap_ill_5um_Grade2]=64G [stratamap_ill_5um_Grade3]=96G )
for TECH in "${SAMPLES[@]}"; do
  # Xenium MM's panel_genes proxy is only 96 genes, below the script's --min-genes default of 200 -
  # override so its shared-panel-genes row (figs 5) isn't silently dropped (see STAGE 4b).
  MIN_GENES=200; [ "$TECH" = "xenium_biomarkers" ] && MIN_GENES=50
  dispatch "extract_counts_$TECH" "${EXTRACT_MEM[$TECH]}" 01:00:00 4 \
    "$R_CMD $CODE_R/extract_count_distributions.R --technology $TECH --min-genes $MIN_GENES --out-dir $REPRO_RESULTS/count_distributions"
done
wait_all

# ============================================================ STAGE 10: N genes / N cells barplot
# data. Feeds fig 5's --ngenes-csv/--ncells-csv.
echo "=== STAGE 10: ngenes/ncells barplot data ==="
dispatch ngenes_ncells 16G 00:15:00 2 \
  "$R_CMD $CODE_R/ngenes_ncells_barplot.R \
    --extract-dir $REPRO_RESULTS/count_distributions \
    --ari-csv $REPRO_RESULTS/separation_ari/separation_ari_comparison.csv \
    --out-dir $REPRO_RESULTS/ngenes_ncells_barplot"
wait_all

# ============================================================ STAGE 11: per-gene-set sensitivity.
# Feeds fig 5's --overall-csv.
echo "=== STAGE 11: sensitivity by gene set ==="
# Exactly the gene sets figure 5 (summary_grid_shared_panel.R) reads for each technology: everyone
# needs "all" (self); atera/xenium_breast100/visiumhd/stratamap x3 have a real shared_panel_genes
# axis; xenium_biomarkers instead needs its own panel_genes (used as its shared-panel proxy);
# xenium_5k needs neither proxy (dropped from panels 3-5 - see summary_grid_shared_panel.R's header).
# --min-count-in-subset 10 matches STAGE 4b/4c's embeddings exactly (same >=10 nCount floor as the
# self/"all" panel's own QC threshold) - so panels 2/3/4/5's shared-panel bars all describe the same
# cell population, not a filtered one for n_cells/ARI and an unfiltered one for nCount/nCount-per-nGenes.
declare -A SENS_GENESETS=( [atera]=all,shared_panel_genes [xenium_breast100]=all,shared_panel_genes \
  [visiumhd]=all,shared_panel_genes [xenium_biomarkers]=all,panel_genes [xenium_5k]=all )
declare -A SENS_MEM=( [atera]=64G [xenium_biomarkers]=64G [xenium_breast100]=64G [xenium_5k]=64G [visiumhd]=64G )
for TECH in "${TECHS[@]}"; do
  dispatch "sensitivity_$TECH" "${SENS_MEM[$TECH]}" 01:30:00 4 \
    "$R_CMD $CODE_R/sensitivity_by_geneset_celltype.R --technology $TECH --sample breast_cancer \
      --qc-rds $(qc_rds_for $TECH) --gene-sets ${SENS_GENESETS[$TECH]} --min-count-in-subset 10 \
      --out-dir $REPRO_RESULTS/sensitivity/$TECH"
done
declare -A SM_SENS_MEM=( [Grade1]=64G [Grade2]=96G [Grade3]=150G )
for GRADE in Grade1 Grade2 Grade3; do
  TECH="stratamap_ill_5um_$GRADE"
  dispatch "sensitivity_sm_$GRADE" "${SM_SENS_MEM[$GRADE]}" 02:00:00 4 \
    "$R_CMD $CODE_R/sensitivity_by_geneset_celltype.R --technology $TECH --sample $GRADE \
      --qc-rds $(qc_rds_for $TECH) --gene-sets all,shared_panel_genes --min-count-in-subset 10 \
      --out-dir $REPRO_RESULTS/sensitivity/$TECH"
done
wait_all
dispatch sensitivity_synthesis 16G 00:15:00 2 \
  "$R_CMD $CODE_R/sensitivity_synthesis.R --sensitivity-dir $REPRO_RESULTS/sensitivity --out-dir $REPRO_RESULTS/sensitivity"
wait_all

# ============================================================ STAGE 12: per-gene sensitivity.
# Feeds fig 8 (--stats-csv). Requires STAGE 5 (spatial metadata) to already exist.
echo "=== STAGE 12: per-gene sensitivity ==="
declare -A GS_MEM=( [atera]=64G [xenium_biomarkers]=64G [xenium_breast100]=64G [xenium_5k]=64G [visiumhd]=96G \
  [stratamap_ill_5um_Grade1]=64G [stratamap_ill_5um_Grade2]=96G [stratamap_ill_5um_Grade3]=150G )
GS_OUT="$REPRO_RESULTS/gene_sensitivity"
mkdir -p "$GS_OUT"
for TECH in "${SAMPLES[@]}"; do
  dispatch "gene_sensitivity_$TECH" "${GS_MEM[$TECH]}" 01:00:00 4 \
    "$R_CMD $CODE_R/gene_sensitivity.R --technology $TECH --out-dir $GS_OUT"
done
wait_all
COMBINED="$GS_OUT/gene_sensitivity_combined.csv"
FIRST=1
for f in "$GS_OUT"/*_gene_sensitivity.csv; do
  if [ "$FIRST" -eq 1 ]; then cat "$f" > "$COMBINED"; FIRST=0; else tail -n +2 "$f" >> "$COMBINED"; fi
done
dispatch gene_sensitivity_synthesis 16G 00:15:00 2 \
  "$R_CMD $CODE_R/gene_sensitivity_synthesis.R --sensitivity-dir $GS_OUT --out-dir $GS_OUT"
wait_all

echo "=== Preprocessing done. Run 03_generate_figures.sh next. ==="
