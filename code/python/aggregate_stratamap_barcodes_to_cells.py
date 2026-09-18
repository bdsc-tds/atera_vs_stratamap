"""
Aggregate StrataMap raw spatial-barcode counts into segmented cells.

StrataMap's raw output (Raw_Matrix_Files/) is a gene x spatial-barcode matrix at sub-cellular
resolution (barcode IDs "SBC:<Y_nm>:<X_nm>", per Illumina's own description: "raw transcript
matrix files (SBC: Y:X in nanometers)"). Segmented cells (from Illumina's ML segmentation) are
provided separately as polygon vertex lists in Segmentations/*_contour_coords_local.csv, in
microns, in the same coordinate frame (no rotation/registration needed - just nm -> um, /1000).

This script rasterizes the cell polygons into a label mask (1 pixel = 1 micron), looks up each raw
barcode's assigned cell via the mask, and sums raw counts per cell via sparse matrix multiplication
to produce a cell x gene matrix, written out in 10x-style mtx format (readable by Seurat::Read10X).

Requires the `pigz` binary on PATH for fast parallel gzip of the (potentially very large) output
matrix; falls back to plain gzip if not found.

Usage:
    python aggregate_stratamap_barcodes_to_cells.py \
        --raw-dir /path/to/<sample>/Raw_Matrix_Files \
        --segmentation-csv /path/to/<sample>/Segmentations/registered_Expanded_5um_cell_contour_coords_local.csv \
        --out-dir /path/to/output/<sample>/raw_aggregated
"""
import argparse
import gzip
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
from PIL import Image, ImageDraw


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_barcodes_um(barcodes_path):
    # Barcode IDs are "SBC:<Y_nm>:<X_nm>" per Illumina docs (Y before X), in nanometers.
    bc = pd.read_csv(barcodes_path, sep=":", header=None, names=["prefix", "y_nm", "x_nm"])
    x_um = bc["x_nm"].to_numpy(dtype=np.float64) / 1000.0
    y_um = bc["y_nm"].to_numpy(dtype=np.float64) / 1000.0
    return x_um, y_um


def load_segmentation(seg_csv):
    seg = pd.read_csv(seg_csv)
    polygons = {
        cell_id: list(zip(g["vertex_x"], g["vertex_y"]))
        for cell_id, g in seg.groupby("cell_id")[["vertex_x", "vertex_y"]]
    }
    return polygons


def rasterize_cells(polygons, xmin, ymin, width, height):
    label_img = Image.new("I", (width, height), 0)
    draw = ImageDraw.Draw(label_img)
    for cell_id, verts in polygons.items():
        px = [(x - xmin, y - ymin) for x, y in verts]
        draw.polygon(px, fill=int(cell_id), outline=int(cell_id))
    return np.asarray(label_img, dtype=np.int64)


def gzip_file(path):
    if shutil.which("pigz"):
        subprocess.run(["pigz", "-f", "-p", "8", str(path)], check=True)
    else:
        with open(path, "rb") as f_in, gzip.open(f"{path}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        Path(path).unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True, help="Raw_Matrix_Files dir (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)")
    ap.add_argument("--segmentation-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--margin-um", type=float, default=10.0, help="Bounding-box margin around segmentation extent")
    args = ap.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log("Loading segmentation polygons")
    polygons = load_segmentation(args.segmentation_csv)
    all_verts = np.array([v for poly in polygons.values() for v in poly])
    seg_xmin, seg_ymin = all_verts[:, 0].min(), all_verts[:, 1].min()
    seg_xmax, seg_ymax = all_verts[:, 0].max(), all_verts[:, 1].max()
    n_cells = len(polygons)
    log(f"{n_cells} segmented cells, extent x=[{seg_xmin:.1f},{seg_xmax:.1f}] y=[{seg_ymin:.1f},{seg_ymax:.1f}] um")

    xmin = seg_xmin - args.margin_um
    ymin = seg_ymin - args.margin_um
    width = int(np.ceil(seg_xmax - xmin + args.margin_um)) + 1
    height = int(np.ceil(seg_ymax - ymin + args.margin_um)) + 1
    log(f"Rasterizing label mask: {width} x {height} px (1um/px)")
    label_arr = rasterize_cells(polygons, xmin, ymin, width, height)

    log("Loading barcode coordinates")
    x_um, y_um = load_barcodes_um(raw_dir / "barcodes.tsv.gz")
    n_barcodes = len(x_um)

    log("Assigning barcodes to cells via label mask lookup")
    px = np.round(x_um - xmin).astype(np.int64)
    py = np.round(y_um - ymin).astype(np.int64)
    in_bounds = (px >= 0) & (px < width) & (py >= 0) & (py < height)
    assigned_cell = np.zeros(n_barcodes, dtype=np.int64)
    assigned_cell[in_bounds] = label_arr[py[in_bounds], px[in_bounds]]
    valid = assigned_cell > 0
    n_valid = valid.sum()
    log(f"{n_valid}/{n_barcodes} barcodes ({100 * n_valid / n_barcodes:.1f}%) fall within a segmented cell")

    cell_ids_sorted = np.sort(np.array(list(polygons.keys())))
    cell_id_to_col_arr = np.full(int(cell_ids_sorted.max()) + 1, -1, dtype=np.int32)
    cell_id_to_col_arr[cell_ids_sorted] = np.arange(n_cells, dtype=np.int32)
    assigned_col = np.full(n_barcodes, -1, dtype=np.int32)
    assigned_col[valid] = cell_id_to_col_arr[assigned_cell[valid]]
    group_col = assigned_col[valid]

    log("Streaming raw gene x barcode matrix, filtering + aggregating per chunk (memory-bounded - "
        "a naive full-matrix read can need >100GB on the largest StrataMap samples)")
    with gzip.open(raw_dir / "matrix.mtx.gz", "rt") as f:
        f.readline()  # %%MatrixMarket...
        f.readline()  # %
        n_genes, _, total_nnz = (int(x) for x in f.readline().split())
    log(f"Raw matrix declared: {n_genes} genes x {n_barcodes} barcodes, nnz={total_nnz}")

    gene_chunks, cell_chunks, count_chunks = [], [], []
    kept = 0
    reader = pd.read_csv(
        raw_dir / "matrix.mtx.gz", sep=" ", skiprows=3, header=None,
        names=["gene_idx", "barcode_idx", "count"],
        dtype={"gene_idx": np.int32, "barcode_idx": np.int32, "count": np.int32},
        chunksize=50_000_000,
    )
    for i, chunk in enumerate(reader):
        cols = assigned_col[chunk["barcode_idx"].to_numpy() - 1]
        mask = cols >= 0
        n_kept = int(mask.sum())
        kept += n_kept
        if n_kept:
            gene_chunks.append(chunk["gene_idx"].to_numpy()[mask] - 1)
            cell_chunks.append(cols[mask])
            count_chunks.append(chunk["count"].to_numpy()[mask])
        log(f"  chunk {i}: {len(chunk)} rows read, {n_kept} kept (cumulative {kept})")

    log("Concatenating filtered entries and building sparse matrix")
    genes_arr = np.concatenate(gene_chunks) if gene_chunks else np.array([], dtype=np.int32)
    cells_arr = np.concatenate(cell_chunks) if cell_chunks else np.array([], dtype=np.int32)
    counts_arr = np.concatenate(count_chunks) if count_chunks else np.array([], dtype=np.int32)
    del gene_chunks, cell_chunks, count_chunks

    cell_gene_counts = scipy.sparse.coo_matrix(
        (counts_arr, (genes_arr, cells_arr)), shape=(n_genes, n_cells)
    ).tocsr()
    cell_gene_counts.sum_duplicates()
    log(f"Aggregated matrix: {cell_gene_counts.shape}, nnz={cell_gene_counts.nnz}")

    log("Computing cell centroids and raw-barcode counts per cell")
    n_raw_per_cell = np.bincount(group_col, minlength=n_cells)
    centroids = {cid: np.mean(polygons[cid], axis=0) for cid in cell_ids_sorted}
    centroid_arr = np.array([centroids[c] for c in cell_ids_sorted])

    log("Writing outputs")
    # dash, not underscore: Seurat::CreateSeuratObject silently replaces "_" with "-" in cell
    # names, which would otherwise desync these barcodes from cell_metadata.parquet downstream.
    barcode_names = [f"cell-{c}" for c in cell_ids_sorted]
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(barcode_names) + "\n")

    features = pd.read_csv(raw_dir / "features.tsv.gz", sep="\t", header=None)
    features.to_csv(out_dir / "features.tsv.gz", sep="\t", header=False, index=False, compression="gzip")

    mtx_path = out_dir / "matrix.mtx"
    scipy.io.mmwrite(str(mtx_path), cell_gene_counts.tocoo(), field="integer")
    gzip_file(mtx_path)

    meta = pd.DataFrame({
        "cell_id": cell_ids_sorted,
        "barcode": barcode_names,
        "x_centroid": centroid_arr[:, 0],
        "y_centroid": centroid_arr[:, 1],
        "n_raw_barcodes": n_raw_per_cell,
    })
    meta.to_parquet(out_dir / "cell_metadata.parquet", index=False)

    log(f"Done. Wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
