"""
Builds the Chromium Flex scRNA-seq reference (janesick2023) used as the RCTD reference for every
technology in this comparison, from two public sources:

  1. GEO GSM7782698 - raw Chromium Flex count matrix (10x-format h5, ~709k barcodes incl. empty
     droplets), from Janesick et al. 2023, Nature Communications 14:8353 ("High resolution mapping
     of the tumor microenvironment using integrated single-cell, spatial and in situ analysis").
  2. Zenodo record 10076046 - per-barcode cell-type annotations for the subset of barcodes that
     passed the paper's own QC/clustering (sheet "scFFPE-Seq" of Cell_Barcode_Type_Matrices.xlsx:
     columns Barcode, Annotation).

The paper's Annotation column (e.g. "Macrophages 1", "DCIS 1", "CD8+ T Cells") is finer-grained
than the Level1/Level2 taxonomy this pipeline's RCTD step expects (Level2 = annotation label,
Level1 = broader class used as RCTD's class_df). reference/janesick_celltype_taxonomy.csv (shipped
in this package) maps every Annotation value to its Level1/Level2 - extracted directly from this
lab's already-built reference object, i.e. this reproduces the exact same reference, not an
approximation.

Cells with no Zenodo annotation (present in the raw GEO matrix but not in the paper's final
annotated set - about 60% of raw barcodes, mostly empty droplets / QC-failed cells) are dropped:
RCTD needs every reference cell to have a real cell type.

Usage: python build_reference.py --out reference/chromium_flex_reference.h5ad
       [--geo-h5 <path>] [--celltype-xlsx <path>] [--taxonomy-csv <path>] [--skip-download]
"""
import argparse
import urllib.request
from pathlib import Path

import pandas as pd
import scanpy as sc

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7782698&format=file&file=GSM7782698_count_raw_feature_bc_matrix.h5"
ZENODO_URL = "https://zenodo.org/records/10076046/files/Cell_Barcode_Type_Matrices.xlsx"


def download(url, dest):
    dest = Path(dest)
    if dest.exists():
        print(f"Already downloaded: {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output .h5ad path")
    ap.add_argument("--geo-h5", default="raw/reference/GSM7782698_count_raw_feature_bc_matrix.h5")
    ap.add_argument("--celltype-xlsx", default="raw/reference/Cell_Barcode_Type_Matrices.xlsx")
    ap.add_argument("--taxonomy-csv", default=str(Path(__file__).resolve().parents[2] / "reference" / "janesick_celltype_taxonomy.csv"),
                     help="Annotation -> Level1/Level2 mapping, ships in this package's reference/ dir")
    ap.add_argument("--skip-download", action="store_true", help="assume --geo-h5/--celltype-xlsx already exist")
    args = ap.parse_args()

    if not args.skip_download:
        download(GEO_URL, args.geo_h5)
        download(ZENODO_URL, args.celltype_xlsx)

    print("Loading raw Chromium Flex counts (GEO GSM7782698)")
    adata = sc.read_10x_h5(args.geo_h5)
    adata.var_names_make_unique()
    print(f"Raw: {adata.shape[0]} barcodes x {adata.shape[1]} genes")

    print("Loading cell-type annotations (Zenodo scFFPE-Seq sheet)")
    annot = pd.read_excel(args.celltype_xlsx, sheet_name="scFFPE-Seq")
    annot = annot.rename(columns={"Barcode": "cell_id", "Annotation": "annotation"}).set_index("cell_id")

    taxonomy = pd.read_csv(args.taxonomy_csv).set_index("annotation")
    annot = annot.join(taxonomy, on="annotation")
    n_unmapped = annot["Level1"].isna().sum() - annot["annotation"].isna().sum()
    if n_unmapped:
        raise SystemExit(f"{n_unmapped} annotation values have no entry in {args.taxonomy_csv} - "
                          "the Zenodo file's Annotation categories may have changed; update the taxonomy CSV")

    common = adata.obs_names.intersection(annot.index)
    print(f"{len(common)} / {adata.shape[0]} barcodes have a cell-type annotation - keeping those only")
    adata = adata[common].copy()
    adata.obs = annot.loc[common, ["annotation", "Level1", "Level2"]].rename(columns={"annotation": "cell_type"})

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"Saved: {args.out} ({adata.shape[0]} cells x {adata.shape[1]} genes)")
    print(adata.obs.groupby(["Level1", "Level2"]).size())


if __name__ == "__main__":
    main()
