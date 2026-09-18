"""
RCTD doublet-mode cell-type annotation for a QCed sample, using the public rctd-py package
(https://github.com/p-gueguen/rctd-py, GPU-accelerated PyTorch reimplementation of the R `spacexr`
RCTD algorithm). Same config for every technology in this comparison: janesick2023 Chromium Flex
reference (see 01_build_reference.py), Level2 annotation, UMI_min=10/counts_MIN=10, class-aware
(class_df from the reference's Level1 column). Output file layout (cell_ids.parquet,
weights.parquet, etc.) matches what code/R/reconstruct_rctd_from_parquet.R expects.

Usage:
    python run_rctd_doublet.py --h5 <sample>_qc_counts.h5 \
        --reference-h5ad <ref>.h5ad --out-dir <rctd_results_dir>
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from rctd import Reference, RCTDConfig, run_rctd

ANNOT_LEVEL = "Level2"
CLASS_LEVEL = "Level1"


def load_reference(reference_h5ad):
    adata_ref = sc.read_h5ad(reference_h5ad)
    adata_ref.var_names_make_unique()
    keep = ~adata_ref.obs[CLASS_LEVEL].isin(["Hybrid / ambiguous"])
    adata_ref = adata_ref[keep].copy()
    return adata_ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", required=True, help="QCed 10x-style h5 counts file")
    ap.add_argument("--reference-h5ad", required=True, help="Chromium reference AnnData (see build_reference.py) - obs must have Level1/Level2 columns")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    save_dir = Path(args.out_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    print("Loading query data")
    query = sc.read_10x_h5(args.h5)
    query.var_names_make_unique()

    print("Loading Chromium reference")
    adata_ref = load_reference(args.reference_h5ad)

    common_genes = adata_ref.var_names.intersection(query.var_names)
    print(f"{len(common_genes)} common genes")

    sc_reference = Reference(adata_ref[:, common_genes].copy(), cell_type_col=ANNOT_LEVEL)
    class_df = dict(zip(adata_ref.obs[ANNOT_LEVEL], adata_ref.obs[CLASS_LEVEL]))

    print("Running RCTD (doublet mode)")
    results = run_rctd(
        query,
        sc_reference,
        mode="doublet",
        config=RCTDConfig(UMI_min=10, counts_MIN=10, class_df=class_df),
    )

    cell_type_names = sc_reference.cell_type_names
    cell_ids = query.obs_names[results.pixel_mask]
    assert len(cell_ids) == results.weights.shape[0]
    assert len(cell_ids) == len(results.spot_class)

    print("Saving results to", save_dir)
    pd.DataFrame({"cell_id": cell_ids}).to_parquet(save_dir / "cell_ids.parquet", index=False)

    pd.DataFrame(results.weights, columns=cell_type_names).to_parquet(save_dir / "weights.parquet")
    pd.DataFrame(results.weights_doublet, columns=["w_1", "w_2"]).to_parquet(save_dir / "weights_doublet.parquet")

    df = pd.DataFrame({
        "cell_id": cell_ids,
        "spot_class": results.spot_class,
        "first_type": results.first_type,
        "second_type": results.second_type,
        "first_class": results.first_class,
        "second_class": results.second_class,
        "min_score": results.min_score,
        "singlet_score": results.singlet_score,
        "first_type_name": [cell_type_names[i] for i in results.first_type],
        "second_type_name": [cell_type_names[i] for i in results.second_type],
    })
    df.to_parquet(save_dir / "spot_results.parquet")

    pd.DataFrame({"cell_id": query.obs_names, "pixel_mask": results.pixel_mask}).to_parquet(save_dir / "pixel_mask.parquet")
    pd.DataFrame({"cell_type_names": cell_type_names}).to_parquet(save_dir / "metadata.parquet")

    with h5py.File(save_dir / "reference_profiles.h5", "w") as f:
        f.create_dataset("profiles", data=sc_reference.profiles)
        f.create_dataset("cell_type_names", data=np.asarray(sc_reference.cell_type_names, dtype="S"))
        f.create_dataset("gene_names", data=np.asarray(sc_reference.gene_names, dtype="S"))

    pd.DataFrame({"DE_genes": sc_reference.get_de_genes()}).to_parquet(save_dir / "de_genes.parquet")

    print("Done:", save_dir.resolve())


if __name__ == "__main__":
    main()
