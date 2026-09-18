"""Build a cell_metadata.parquet for VisiumHD's cell level from Space Ranger 4.1.0's
segmented_outputs bundle. The bundle already ships a ready-made cell x gene matrix
(filtered_feature_cell_matrix.h5/.mtx.gz - standard 10x format, no conversion needed for counts) but
no metadata parquet the way Xenium-family cells.parquet does - geometry only exists as two separate
GeoJSON files (cell_segmentations.geojson, nucleus_segmentations.geojson), joined by an integer
cell_id that maps to the matrix barcodes as "cellid_{cell_id:09d}-1".

The older binned_outputs.tar.gz bundle (square_002um/008um/016um) is a separate, bin-only download;
segmented_outputs is 10x's own newer built-in cell-segmentation pipeline (Space Ranger >=4.1.0).

IMPORTANT: cell_segmentations.geojson's
coordinates are in FULL-RESOLUTION IMAGE PIXELS, not microns - confirmed by converting the observed
x-span (~40,416 px) via spatial/scalefactors_json.json's microns_per_pixel (0.2737094905827413) to
~11,062um, matching the "11mm" capture area almost exactly. The original version of this script used
the raw GeoJSON values directly as if they were already microns, inflating x_centroid/y_centroid by
~1/0.2737~=3.65x (linear) and cell_area/nucleus_area by ~13.3x (quadratic, since area scales as the
square of the linear factor) - this fed wrong-scale coordinates into the spatial spillover network
(SPLIT::build_spatial_network's rad_pruning is a real-micron radius) and made "VisiumHD cells are ~10x
bigger" look like a real finding when it was this unit bug. MICRONS_PER_PIXEL below is read from that
scale factor file at runtime, not hardcoded, so it's tied to this dataset's own calibration.

Usage: python3 convert_visiumhd_segmented_cells.py <segmented_outputs_dir> <output_dir> <scalefactors_json_path>
"""
import sys
import json
import numpy as np
import pandas as pd


def shoelace_area(coords):
    """coords: list of [x, y] pairs, first ring only (exterior) - all polygons here are simple,
    no holes, per direct inspection of the GeoJSON."""
    x = np.array([p[0] for p in coords])
    y = np.array([p[1] for p in coords])
    return 0.5 * abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def load_geojson_polygons(path, area_col):
    with open(path) as f:
        gj = json.load(f)
    rows = []
    for feat in gj["features"]:
        props = feat["properties"]
        coords = feat["geometry"]["coordinates"][0]  # exterior ring
        rows.append({"cell_id_int": props["cell_id"], area_col: shoelace_area(coords)})
    return pd.DataFrame(rows)


def main():
    seg_dir, out_dir, scalefactors_path = sys.argv[1], sys.argv[2], sys.argv[3]
    import os
    os.makedirs(out_dir, exist_ok=True)

    with open(scalefactors_path) as f:
        microns_per_pixel = json.load(f)["microns_per_pixel"]
    print(f"microns_per_pixel = {microns_per_pixel} (from {scalefactors_path})")

    print("Parsing cell_segmentations.geojson ...")
    with open(f"{seg_dir}/cell_segmentations.geojson") as f:
        cell_gj = json.load(f)
    cell_rows = []
    for feat in cell_gj["features"]:
        props = feat["properties"]
        coords = feat["geometry"]["coordinates"][0]
        cell_rows.append({
            "cell_id_int": props["cell_id"],
            "x_centroid": props["cell_centroid"][0],
            "y_centroid": props["cell_centroid"][1],
            "nucleus_x_centroid": props["nucleus_centroid"][0],
            "nucleus_y_centroid": props["nucleus_centroid"][1],
            "cell_area": shoelace_area(coords),
        })
    cell_df = pd.DataFrame(cell_rows)
    print(f"{len(cell_df)} cells")

    print("Parsing nucleus_segmentations.geojson ...")
    nuc_df = load_geojson_polygons(f"{seg_dir}/nucleus_segmentations.geojson", "nucleus_area")
    print(f"{len(nuc_df)} nuclei")

    merged = cell_df.merge(nuc_df, on="cell_id_int", how="left")
    merged["cell_id"] = merged["cell_id_int"].apply(lambda i: f"cellid_{i:09d}-1")
    merged = merged.drop(columns=["cell_id_int"])

    # Pixels -> microns: linear columns scale by microns_per_pixel, area columns by its square.
    for col in ["x_centroid", "y_centroid", "nucleus_x_centroid", "nucleus_y_centroid"]:
        merged[col] = merged[col] * microns_per_pixel
    for col in ["cell_area", "nucleus_area"]:
        merged[col] = merged[col] * (microns_per_pixel ** 2)

    merged = merged[["cell_id", "x_centroid", "y_centroid", "cell_area",
                      "nucleus_x_centroid", "nucleus_y_centroid", "nucleus_area"]]

    out_path = f"{out_dir}/cell_metadata.parquet"
    merged.to_parquet(out_path, index=False)
    print(f"Saved {out_path}: {len(merged)} cells")
    print(merged.describe())


if __name__ == "__main__":
    main()
