"""
prepare_wfr_data.py
--------------------
Combines the four Wildfire Risk MCDA GeoPackage outputs (Final WFR,
Wildfire Hazard, Wildfire Exposure, Social Vulnerability Index — each
a zonal-statistics join onto the 55 San Diego County Populated Places)
into a single GeoJSON for mapping, with each submodel's mean/rank
carried as prefixed properties.

This version uses gpkg_reader.py (pure Python + a hand-rolled Albers
Conic inverse projection) so it runs without GDAL/geopandas installed.
If you have GeoPandas + pyproj available, this is equivalent to:

    gdf = gpd.read_file(path).to_crs(4326)

...for each layer, followed by a merge on GEOID.

Usage:
    python prepare_wfr_data.py --uploads /path/to/gpkg/files \
        --output ../data/wfr_populated_places.geojson
"""

import argparse
import json
import os

from gpkg_reader import read_gpkg_layer

LAYERS = {
    "wfr": "Final_WFR_Populated_Places_25_05_09_11_12.gpkg",
    "wfh": "WFH_Populated_Places_25_05_09_11_12.gpkg",
    "wfe": "WFE_Populated_Places_25_05_09_11_12.gpkg",
    "svi": "SVI_Populated_Places_25_05_09_11_12.gpkg",
}


def multipolygon_to_geojson_coords(multipoly):
    # multipoly: list[polygon], polygon: list[ring], ring: list[(lon,lat)]
    return [[[[lon, lat] for lon, lat in ring] for ring in poly] for poly in multipoly]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uploads", default="../raw")
    ap.add_argument("--output", default="../data/wfr_populated_places.geojson")
    args = ap.parse_args()

    layer_data = {}
    base_geom_by_geoid = {}
    for key, fname in LAYERS.items():
        path = os.path.join(args.uploads, fname)
        feats = read_gpkg_layer(path)
        by_geoid = {}
        for attrs, geom in feats:
            geoid = attrs["GEOID"]
            by_geoid[geoid] = attrs
            if key == "wfr":
                base_geom_by_geoid[geoid] = (attrs["NAME"], geom)
        layer_data[key] = by_geoid

    features = []
    for geoid, (name, geom) in base_geom_by_geoid.items():
        props = {"GEOID": geoid, "NAME": name}
        for key in LAYERS:
            attrs = layer_data[key].get(geoid, {})
            props[f"{key}_mean"] = round(attrs.get("mean", 0), 4)
            props[f"{key}_rank"] = attrs.get("rank")
            props[f"{key}_max"] = round(attrs.get("max", 0), 4) if attrs.get("max") is not None else None
            props[f"{key}_std"] = round(attrs.get("std", 0), 4) if attrs.get("std") is not None else None

        features.append({
            "type": "Feature",
            "properties": props,
            "geometry": {
                "type": "MultiPolygon",
                "coordinates": multipolygon_to_geojson_coords(geom),
            },
        })

    fc = {"type": "FeatureCollection", "features": features}
    with open(args.output, "w") as f:
        json.dump(fc, f)

    print(f"Wrote {len(features)} populated-place features to {args.output}")


if __name__ == "__main__":
    main()
