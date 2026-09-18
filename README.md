# Reproducing the Atera vs Stratamap analysis posted at https://lnkd.in/p/gHZ84RsX

A self-contained package to reproduce the 8 figures from that post: a spatial-transcriptomics
technology comparison (Atera, VisiumHD, Xenium MM/biomarkers panel, Xenium v1/breast+100 addon
panel, Xenium 5k, StrataMap - 6 technologies, 8 cell-level samples, all annotated against the same
10x Chromium Flex scRNA-seq reference). Nothing in this package points at any specific machine, lab,
or institution: raw data comes from public download links or a documented manual step (see
`00_download_raw_data.md`), and every other path is an environment variable you set (see `code/lib.sh`).

**[Read the full analysis report](04_analysis_report.md)** - narrates all 8 figures in order, with
the exact code to reproduce each one.

## The 8 figures

1. One spatial scatter per sample, colored by cell type
2. Each technology's own UMAP
3. Cell-type composition per sample
4. % cells retained vs. nCount QC floor, by cell type
5. N genes / N cells / nCount / ARI, self vs. shared-panel-genes
6. RCTD secondary-type-weight distribution per sample
7. w2-vs-spatial-neighborhood cosine similarity, ranked
8. Per-gene sensitivity, own panel vs. shared gene axis

Exact output filenames and reproduction commands for each are in `04_analysis_report.Rmd` and
`03_generate_figures.sh`.

## Layout

- **`00_download_raw_data.md`** - every raw dataset + where to get it. All public except
  StrataMap (Illumina BaseSpace, browser download).
- **`environment.yml`** - the R + Python environment (conda). Two R packages aren't on
  CRAN/Bioconductor and need a separate install step - see the comment at the top of
  `02_preprocess.sh`.
- **`02_preprocess.sh`** - raw data -> every processed artifact the figures need: per-technology
  ingestion (including StrataMap's raw-barcode-to-cell aggregation) -> QC -> RCTD cell-type
  annotation -> gene-set embeddings -> spatial metadata -> spillover -> composition/sensitivity/ARI
  summaries.
- **`03_generate_figures.sh`** - the 8 exact plotting commands, run after `02_preprocess.sh`.
- **`04_analysis_report.Rmd`** - source for the report; **`04_analysis_report.md`** is the rendered,
  GitHub-browsable version (images in `report_figures/`), **`04_analysis_report.html`** a
  self-contained standalone version (same content, everything embedded in one file). Re-knit either
  with `rmarkdown::render("04_analysis_report.Rmd", output_format = "github_document")` or
  `"html_document"`. Degrades gracefully: shows the real rendered PNG if `03_generate_figures.sh`
  has run, otherwise rebuilds a plain version from whatever CSVs `02_preprocess.sh` has already
  produced (figure 2's UMAPs are the one exception - no CSV fallback is possible for those).
- **`code/R/`, `code/python/`** - every script the two orchestration scripts above call.
- **`reference/`** - small, non-sensitive lookup tables the pipeline needs: gene panel lists per
  technology, a gene biotype (protein-coding/lncRNA/etc.) table, the shared cell-type color
  palette, and the Chromium reference's cell-type taxonomy (see below).

## How to run

```bash
conda env create -f environment.yml
conda activate spatial-tech-comparison
R -e 'remotes::install_github("dmcable/spacexr")'   # RCTD (R side, used to reconstruct results for spillover.R)
R -e 'remotes::install_github("bdsc-tds/SPLIT")'    # spatial-neighborhood spillover methodology

# 1. Get the raw data - see 00_download_raw_data.md. Save everything under one directory, e.g.:
export REPRO_DATA=/path/with/room/for/large/data
mkdir -p "$REPRO_DATA/raw"
# ... download each dataset into $REPRO_DATA/raw/<technology>/ as documented ...

# 2. Preprocess (raw -> every intermediate artifact). This is the expensive part - RCTD and the
#    gene-set embeddings need real memory (32-250GB depending on technology/sample size, StrataMap
#    Grade3's ~700k cells being the largest) - see the SLURM notes below if you have cluster access.
bash 02_preprocess.sh

# 3. Generate the figures (cheap - seconds to minutes).
bash 03_generate_figures.sh
```

Figures land in `$REPRO_RESULTS/figures/` (defaults to `./results/figures` next to this README).

### Running on a SLURM cluster

Set `USE_SLURM=true` and `SLURM_ACCOUNT=<your account>` before running `02_preprocess.sh` - every
compute-heavy step then runs as a separate `sbatch` job instead of directly on whatever machine
you launched the script from. Without SLURM, `02_preprocess.sh` runs everything **sequentially** by
design: several of these steps need 32-250GB of RAM, and running multiple at once on a single
workstation can exhaust its memory even though any one of them alone is fine.

## Notes / scope

- **StrataMap's fresh-frozen vs. everyone-else's FFPE tissue prep is a real confound baked into
  every one of these 8 figures**
- **RCTD** runs via `rctd-py` (https://github.com/p-gueguen/rctd-py), a public GPU-accelerated
  PyTorch reimplementation - not the original lab's private tooling. Every run in this package used
  plain CPU nodes (no `--gres=gpu` anywhere in `02_preprocess.sh`), including VisiumHD's ~850k cells
  in well under an hour - a GPU isn't necessary for datasets this size.
