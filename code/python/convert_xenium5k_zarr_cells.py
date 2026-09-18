"""Convert the Xenium 5k breast dataset's Explorer zarr bundle (cell_feature_matrix.zarr.zip +
cells.zarr.zip) to the same 10x-style mtx + parquet metadata shape the other Xenium-family
technologies (Atera, biomarkers, breast+100) ship natively at cell level - so downstream Stage B/C
scripts can treat every technology's cell level identically.

This dataset ("xe_outs" bundle, Xenium Explorer's zarr export) has no cell_feature_matrix.h5/
cells.parquet. cells.zarr.zip's `cell_summary` array's own `.attrs['column_names']` (not the
group's attrs - easy to miss) documents its 8 columns: cell_centroid_x/y, cell_area,
nucleus_centroid_x/y, nucleus_area, z_level, nucleus_count.

Bin-level (needs the transcripts.zarr.zip grid-pyramid's finest level, concatenated across tiles) and
nucleus-level (needs per-transcript positions mapped through cells.zarr.zip's pixel-level
segmentation mask via `homogeneous_transform`, since there's no per-transcript overlaps_nucleus flag
in this bundle) are NOT handled here - cell level only.

Usage: python3 convert_xenium5k_zarr_cells.py <xe_outs_dir> <output_dir>
"""
import sys
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import zarr


def main():
    xe_outs_dir, out_dir = sys.argv[1], sys.argv[2]
    import os
    os.makedirs(out_dir, exist_ok=True)

    print("Loading cell_feature_matrix.zarr.zip ...")
    cfm = zarr.open(zarr.storage.ZipStore(f"{xe_outs_dir}/cell_feature_matrix.zarr.zip", mode="r"), mode="r")
    cf = cfm["cell_features"]
    attrs = dict(cf.attrs)
    n_features = attrs["number_features"]
    n_cells = attrs["number_cells"]
    feature_keys = np.array(attrs["feature_keys"])   # gene symbols
    feature_ids = np.array(attrs["feature_ids"])     # Ensembl IDs
    feature_types = np.array(attrs["feature_types"])
    print(f"{n_features} features x {n_cells} cells, feature_types: "
          f"{pd.Series(feature_types).value_counts().to_dict()}")

    # csc group is already (features x cells) in CSC form - columns=cells, matches Seurat convention
    csc = cf["csc"]
    mat = sp.csc_matrix(
        (csc["data"][:], csc["indices"][:], csc["indptr"][:]),
        shape=(n_features, n_cells),
    )
    print(f"Reconstructed CSC matrix: {mat.shape}, nnz={mat.nnz}")

    cfm_cell_id = cf["cell_id"][:]  # (n_cells, 2) uint32

    print("Loading cells.zarr.zip ...")
    cz = zarr.open(zarr.storage.ZipStore(f"{xe_outs_dir}/cells.zarr.zip", mode="r"), mode="r")
    cells_cell_id = cz["cell_id"][:]  # (n_cells, 2) uint32
    # Sanity check: cell_feature_matrix's cell axis order must match cells.zarr's cell_id order -
    # assumed positional alignment (both length n_cells, both from the same XOA pipeline run), not
    # independently guaranteed by any join key - verify before trusting it.
    assert np.array_equal(cfm_cell_id, cells_cell_id), (
        "cell_id order mismatch between cell_feature_matrix.zarr and cells.zarr - "
        "positional alignment assumption is WRONG, need an explicit join instead"
    )
    print("Confirmed: cell_feature_matrix.zarr and cells.zarr share identical cell_id order (positional join is safe)")

    cs_attrs = dict(cz["cell_summary"].attrs)
    col_names = cs_attrs["column_names"]
    cell_summary = pd.DataFrame(cz["cell_summary"][:], columns=col_names)
    print("cell_summary columns:", col_names)

    # Synthesize a string cell_id from the (row, suffix) uint32 pair - doesn't need to match 10x's
    # own internal string encoding, just needs to be unique and stable within this dataset.
    cell_id_str = [f"5k-{a}-{b}" for a, b in cells_cell_id]

    # IMPORTANT: feature_type=='aggregate_gene' (symbol "Total transcripts") is a precomputed
    # per-cell SUM of all real gene rows, not a real feature - confirmed by direct comparison
    # (its row-sum across all cells exactly equals the sum of all feature_type=='gene' rows).
    # Summing mat over ALL rows double-counts gene signal via this row (verified: gave 154.6M vs. a
    # correct 71.6M gene-only total, implausibly higher than experiment.xenium's own
    # num_transcripts=92.8M). Restrict total_counts/n_genes_detected to real gene rows only - keep
    # the aggregate/control rows in the matrix itself for downstream flexibility, just don't sum
    # over them here.
    gene_mask = feature_types == "gene"
    total_counts = np.asarray(mat[gene_mask].sum(axis=0)).ravel()
    n_genes_detected = np.asarray((mat[gene_mask] > 0).sum(axis=0)).ravel()

    cell_metadata = pd.DataFrame({
        "cell_id": cell_id_str,
        "x_centroid": cell_summary["cell_centroid_x"],
        "y_centroid": cell_summary["cell_centroid_y"],
        "cell_area": cell_summary["cell_area"],
        "nucleus_x_centroid": cell_summary["nucleus_centroid_x"],
        "nucleus_y_centroid": cell_summary["nucleus_centroid_y"],
        "nucleus_area": cell_summary["nucleus_area"],
        "nucleus_count": cell_summary["nucleus_count"],
        "total_counts": total_counts,
        "n_genes_detected": n_genes_detected,
    })

    print(f"Writing outputs to {out_dir} ...")
    cell_metadata.to_parquet(f"{out_dir}/cell_metadata.parquet", index=False)

    features_df = pd.DataFrame({
        "feature_id": feature_ids,
        "gene_symbol": feature_keys,
        "feature_type": feature_types,
    })
    features_df.to_csv(f"{out_dir}/features.tsv.gz", sep="\t", header=False, index=False,
                        columns=["feature_id", "gene_symbol", "feature_type"])
    with gzip.open(f"{out_dir}/barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(cell_id_str) + "\n")
    # scipy.io.mmwrite on an unsigned-int dtype (uint32, from the zarr source arrays) writes a
    # MatrixMarket "unsigned-integer" field type, which R's Matrix::readMM (used internally by
    # Seurat::Read10X) does not recognize ("element type 'unsigned-integer' not recognized" -
    # confirmed by hitting this directly in scripts/02_QC_seurat.R). Cast to signed int32 first so
    # the header says plain "integer", which readMM does handle. Counts are small (max seen here:
    # low thousands per cell), no overflow risk from the uint32->int32 cast.
    mat = mat.astype(np.int32)
    sio.mmwrite(f"{out_dir}/matrix.mtx", mat)
    import subprocess
    subprocess.run(["gzip", "-f", f"{out_dir}/matrix.mtx"], check=True)

    print("Done.")
    print(f"  cell_metadata.parquet: {len(cell_metadata)} cells")
    print(f"  features.tsv.gz: {len(features_df)} features")
    print(f"  matrix.mtx.gz: {mat.shape}, nnz={mat.nnz}")


if __name__ == "__main__":
    main()
