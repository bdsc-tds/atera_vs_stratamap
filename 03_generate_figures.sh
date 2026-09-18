#!/bin/bash
# Regenerates the 8 specific figures, assuming 02_preprocess.sh has already produced every input
# file each one reads. Every step (including plotting) goes through dispatch() from code/lib.sh,
# same as 02_preprocess.sh - set USE_SLURM=true to run these as SLURM jobs too, rather than on
# whatever node you launched the script from.
#
# Figure 1 was originally produced with hand-tuned --point-size/--alpha for a "bigger dots, more
# visible at a glance" look and a taller-than-wide per-panel aspect ("asp2" = panel-height/
# panel-width = 2). The aspect is exactly reconstructed (verified by back-solving the original
# PNG's pixel dimensions: 3600x900px / 300dpi = 12x3in = 8 panels * 1.5in x 3in), but the exact
# point-size/alpha values used originally weren't preserved anywhere retrievable - the values below
# are a reasonable starting point per spatial_grid.R's own header comment; adjust to taste.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/code/lib.sh"

OUT="$REPRO_RESULTS/figures"

dispatch fig1_spatial_grid 32G 00:20:00 2 \
  "$R_CMD $CODE_R/spatial_grid.R --out-dir $OUT/spatial --format png --dpi 300 \
    --point-size 0.3 --alpha 0.7 --panel-width 1.5 --panel-height 3 --out-name spatial_grid_all_bigdots_asp2"

dispatch fig2_umap_by_geneset 16G 00:20:00 2 \
  "$R_CMD $CODE_R/umap_by_geneset.R --out-dir $OUT/umap"
dispatch fig2_umap_self_grid 16G 00:20:00 2 \
  "$R_CMD $CODE_R/umap_self_grid.R --out-dir $OUT/umap --out-name umap_self_grid_all"

# composition_summary.csv already exists from 02_preprocess.sh's STAGE 7 (celltype_composition.R) -
# no need to recompute it here, composition_grid.R just draws figure 3 from it directly.
dispatch fig3_composition_grid 16G 00:15:00 2 \
  "$R_CMD $CODE_R/composition_grid.R --summary-csv $REPRO_RESULTS/composition/composition_summary.csv --out-dir $OUT/composition/sorted_by_abundance"

dispatch fig4_threshold_sensitivity 16G 00:20:00 2 \
  "$R_CMD $CODE_R/celltype_threshold_sensitivity.R --extract-dir $REPRO_RESULTS/count_distributions --out-dir $OUT/threshold_sensitivity --layout row_tall"

dispatch fig5_summary_grid 16G 00:20:00 2 \
  "$R_CMD $CODE_R/summary_grid_shared_panel.R \
    --ngenes-csv $REPRO_RESULTS/ngenes_ncells_barplot/ngenes_barplot.csv \
    --ncells-csv $REPRO_RESULTS/ngenes_ncells_barplot/ncells_barplot.csv \
    --overall-csv $REPRO_RESULTS/sensitivity/sensitivity_overall_all_technologies.csv \
    --ari-csv $REPRO_RESULTS/separation_ari/separation_ari_comparison.csv \
    --out-dir $OUT/summary_grid"

dispatch fig6_w2_boxplot 16G 00:15:00 2 \
  "$R_CMD $CODE_R/w2_boxplot.R --cache-root $REPRO_RESULTS/spillover_cache --out-dir $OUT/w2"

dispatch fig7_spillover_metrics 16G 00:20:00 2 \
  "$R_CMD $CODE_R/spillover_metrics.R --cache-root $REPRO_RESULTS/spillover_cache --out-dir $OUT/spillover"

dispatch fig8_gene_sensitivity_scatter 16G 00:15:00 2 \
  "$R_CMD $CODE_R/gene_sensitivity_scatter.R --stats-csv $REPRO_RESULTS/gene_sensitivity/csv/gene_sensitivity_stats.csv --out-dir $OUT/gene_sensitivity_scatter"

wait_all

echo "Done. Figures:"
echo "  1) $OUT/spatial/spatial_grid_all_bigdots_asp2.png"
echo "  2) $OUT/umap/umap_self_grid_all.png (+ per-sample $OUT/umap/umap_self_<tech>.png)"
echo "  3) $OUT/composition/sorted_by_abundance/composition_grid_all_pct_linear_fixedy_tall.png"
echo "  4) $OUT/threshold_sensitivity/celltype_threshold_sensitivity_tall.png"
echo "  5) $OUT/summary_grid/summary_grid_shared_panel.png"
echo "  6) $OUT/w2/w2_boxplot_by_score.png"
echo "  7) $OUT/spillover/spillover_cosine.png"
echo "  8) $OUT/gene_sensitivity_scatter/scatter_self_vs_shared_panel_genes_transcripts_per_cell_median_log10.png"
