"""
build_wfr_map.py
------------------
Builds the Wildfire Risk MCDA choropleth webmap (55 San Diego County
Populated Places, 4 submodels) using leafmap + GeoPandas.

Run this in an environment with the geospatial stack installed:
    pip install leafmap geopandas mapclassify

Reads the four GeoPackages directly (no manual reprojection needed —
GeoPandas + pyproj handle the EPSG:3310 -> EPSG:4326 transform), merges
them on GEOID, and renders a submodel-togglable choropleth.

Usage:
    python build_wfr_map.py --uploads /path/to/gpkg/files \
        --output ../maps/wfr_populated_places_map.html
"""

import argparse
import os

import geopandas as gpd
import leafmap.foliumap as leafmap

LAYERS = {
    "wfr": "Final_WFR_Populated_Places_25_05_09_11_12.gpkg",
    "wfh": "WFH_Populated_Places_25_05_09_11_12.gpkg",
    "wfe": "WFE_Populated_Places_25_05_09_11_12.gpkg",
    "svi": "SVI_Populated_Places_25_05_09_11_12.gpkg",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uploads", default=".")
    ap.add_argument("--output", default="../maps/wfr_populated_places_map.html")
    args = ap.parse_args()

    merged = None
    for key, fname in LAYERS.items():
        gdf = gpd.read_file(os.path.join(args.uploads, fname)).to_crs(4326)
        gdf = gdf[["GEOID", "NAME", "mean", "rank", "max", "std", "geometry"]].rename(
            columns={"mean": f"{key}_mean", "rank": f"{key}_rank",
                     "max": f"{key}_max", "std": f"{key}_std"}
        )
        if merged is None:
            merged = gdf
        else:
            merged = merged.merge(
                gdf.drop(columns="geometry"), on=["GEOID", "NAME"], how="left"
            )

    m = leafmap.Map(center=[32.95, -116.85], zoom=9)
    m.add_data(
        merged,
        column="wfr_mean",
        cmap="OrRd",
        legend_title="Final Wildfire Risk (mean)",
        layer_name="Wildfire Risk MCDA",
        info_mode="on_hover",
    )
    # To add hazard / exposure / SoVI as togglable layers instead of a single
    # column, repeat m.add_data(...) with column="wfh_mean" etc., each as its
    # own named layer, then m.add_layer_control().

    m.to_html(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
