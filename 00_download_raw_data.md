# Raw data manifest

Everything below is public except StrataMap. Each 10x Genomics dataset page is behind a JS/bot-check
that blocks scripted fetches (`curl`, browser automation, etc. all get blocked/429'd) - so those
need one manual browser download each. VisiumHD and the Chromium reference are the two exceptions:
their files sit on plain CDNs with no bot-check, so `02_preprocess.sh` and
`code/python/build_reference.py` fetch them automatically.

Save everything under a single directory (`$REPRO_DATA/raw/`, see `02_preprocess.sh`'s config
block) using the subfolder names below.

## 1. Chromium Flex scRNA-seq reference (janesick2023)

Used as the RCTD reference for every technology (Level2 = annotation label, Level1 = class_df).
Built automatically by `code/python/build_reference.py` from two public sources - no manual
download needed:

- Raw counts (GEO): https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782698&format=file&file=GSM7782698_count_raw_feature_bc_matrix.h5
- Cell-type annotations (Zenodo): https://zenodo.org/records/10076046/files/Cell_Barcode_Type_Matrices.xlsx
- Paper: Janesick et al. 2023, *Nature Communications* 14:8353 - "High resolution mapping of the
  tumor microenvironment using integrated single-cell, spatial and in situ analysis"
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC10730913/)
- Companion code: https://github.com/10XGenomics/janesick_nature_comms_2023_companion

The Zenodo file's per-cell "Annotation" column is finer-grained than the Level1/Level2 taxonomy
this pipeline needs; `reference/janesick_celltype_taxonomy.csv` (shipped in this package) maps
every Annotation value to Level1/Level2 - see `build_reference.py`'s header for how that mapping
was derived.

## 2. Atera (Xenium WTA Preview), FFPE human breast cancer

- Dataset page: https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer
- Download the "Xenium Output Bundle" (standard `outs/` layout: `cell_feature_matrix.h5`,
  `cells.parquet`, `cell_boundaries.parquet`, `transcripts.parquet`, ...).
- Save/unzip to: `$REPRO_DATA/raw/atera/outs/`

## 3. Xenium targeted - Biomarkers panel ("Xenium MM")

- Dataset page: https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers
- Download the Xenium Output Bundle.
- Save/unzip to: `$REPRO_DATA/raw/xenium_biomarkers/outs/`

## 4. Xenium targeted - Breast + 100 addon panel ("Xenium v1")

- Dataset page: https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-with-custom-add-on-panel-1-standard
- Download the Xenium Output Bundle.
- Save/unzip to: `$REPRO_DATA/raw/xenium_breast100/outs/`

## 5. Xenium Prime 5K ("Xenium 5k"), FFPE human breast cancer

- Dataset page: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer
- Download the Xenium Output Bundle - **this one ships as a zarr-based Explorer bundle**
  (`cell_feature_matrix.zarr.zip`, `cells.zarr.zip`, `transcripts.zarr.zip`), not the parquet/h5
  `outs/` shape the other 3 Xenium-family datasets use. `02_preprocess.sh` converts it with
  `code/python/convert_xenium5k_zarr_cells.py`.
- Save/unzip to: `$REPRO_DATA/raw/xenium_5k/xe_outs/`
- A lower-res, cropped SpatialData mirror also exists on Zenodo if useful for quick exploration
  (not what this pipeline uses): https://zenodo.org/records/20180722

## 6. VisiumHD, 11mm human breast cancer

- Dataset page (for reference/ToS): https://www.10xgenomics.com/datasets (search "Visium HD 11mm
  Human Breast Cancer")
- Direct CDN base (no bot-check, reachable via `curl -I` -> HTTP 200):
  `https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_11mm_Human_Breast_Cancer/`
- File actually needed by this pipeline (cell level only):
  `Visium_HD_11mm_Human_Breast_Cancer_segmented_outputs.tar.gz` (Space Ranger 4.1.0's built-in cell
  segmentation - ships `filtered_feature_cell_matrix.h5` + segmentation GeoJSONs +
  `spatial/scalefactors_json.json`)
- `02_preprocess.sh` downloads this directly with `curl`.
- Save to: `$REPRO_DATA/raw/visiumhd/`

## 7. StrataMap - human breast cancer demo set (DCIS + IDC grades 1-3)

**Not a plain download - needs an Illumina BaseSpace DataCentral account.**

- Portal: https://basespace.illumina.com/datacentral
- The "authenticated API" route (BaseSpace CLI / OAuth device-flow / app access-token) does not
  work on a trial account - trial accounts have no API access, both the OAuth device-flow and an
  app access-token are rejected. If you have a paid/full BaseSpace account, the API route may work
  - see Illumina's BaseSpace Sequence Hub API docs.
- **What actually works**: sign in to BaseSpace DataCentral, find the StrataMap demo breast cancer
  dataset, and use Download > Analysis (the only enabled option on a trial account) via the
  browser. This gives you, per sample: `Raw_Matrix_Files/{matrix.mtx.gz,features.tsv.gz,
  barcodes.tsv.gz}`, `Images/registered.ome.tiff`, `Segmentations/registered_{Expanded_5um_cell,
  nuclei}_contour_coords_local.csv`.
- 3 samples: `DCIS_IDC_Grade1-77125049`, `IDC_Grade2-77036971`, `IDC_Grade3-77134058`.
- Save to: `$REPRO_DATA/raw/stratamap/<sample>/`
- **Important confound**: StrataMap only supports fresh-frozen (FF) tissue (FFPE support not yet
  available at time of writing) - every other technology here is FFPE. This is baked into the
  comparison, not a side hypothesis.

`02_preprocess.sh` includes the full StrataMap ingestion (raw-barcode-to-cell aggregation, QC,
RCTD) - once you have the BaseSpace download, everything downstream is automated the same as every
other technology.
