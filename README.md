# Wildfire Risk MCDA — San Diego County

A GIS multi-criteria decision analysis (MCDA) of wildfire risk across San Diego County. It combines a **Wildfire Hazard (WFH)**, a **Wildfire Exposure (WFE)** and a **Social Vulnerability Index (SVI)** submodel into a composite **Wildfire Risk (WFR)** surface, then summarizes the results for 55 populated places.

**[View the interactive webmap →](https://mtwyman6243.github.io/research.html)** · Part of [Matt Twyman's GIS portfolio](https://mtwyman6243.github.io)

## About the project

This is the M.S. thesis research (San Diego State University, Department of Geography). It was funded by NASA's Equity and Environmental Justice (EEJ) Program, Award 21-EEJ21-0037. Co-authors: Douglas Stow, Megan Jennings, Piotr Jankowski (SDSU) and Gernot Paulus (Carinthia University of Applied Sciences).

- Submodel criteria were weighted with the Analytic Hierarchy Process (AHP) and aggregated with a Weighted Linear Combination at 30 m resolution.
- Variance-based sensitivity and uncertainty analysis (Monte Carlo simulation, Sobol indices) tested how robust the model is.
- Boulevard and the southeastern Tribal communities came out highest-risk. WUI classification and structure density were the strongest drivers of exposure.

## What's in this repo

```
scripts/
  prepare_wfr_data.py   # merges the 4 submodel GeoPackages -> one GeoJSON (pure Python)
  build_wfr_map.py      # GeoPandas + leafmap version of the webmap build
  gpkg_reader.py        # minimal GeoPackage reader + hand-rolled Albers (EPSG:3310) inverse
  geotiff_reader.py     # minimal GeoTIFF reader (struct + numpy)
  raster_to_web.py      # reprojects 30 m rasters to WGS84 PNG overlays for Leaflet
  austria_projection.py # shared projection helper imported by raster_to_web.py
data/
  wfr_populated_places.geojson   # 55 populated places with all submodel scores/ranks
maps/
  wfr_populated_places_map.html  # standalone Leaflet webmap
  assets/                        # raster overlays (PNG) + WGS84 bounds (JSON)
```

## The pipeline

The preprocessing runs offline once. Its output is static GeoJSON and PNG files, so the map needs no server.

1. **Vectors:** `prepare_wfr_data.py` reads the four zonal-statistics GeoPackages, reprojects them from NAD83 California Albers (EPSG:3310) to WGS84, and merges them on `GEOID`.
2. **Rasters:** `raster_to_web.py` reads each GeoTIFF and resamples it onto a WGS84 grid using backward mapping: each output cell is forward-projected into Albers to look up its source pixel. It then writes a colorized PNG plus its bounds.
3. **Map:** a Leaflet page with a submodel toggle (WFR / WFH / WFE / SVI) and a switch between the place-level choropleth and the continuous raster.

The projection and file-format code is written from scratch (Snyder 1987 formulas, raw GeoPackage WKB and TIFF IFD parsing) rather than using GDAL or pyproj. `build_wfr_map.py` shows the equivalent GeoPandas approach.

## Running it

The raw inputs (GeoPackages and GeoTIFFs) aren't committed. Put them in a `raw/` folder, then:

```bash
pip install -r requirements.txt
cd scripts
python prepare_wfr_data.py --uploads ../raw --output ../data/wfr_populated_places.geojson
python raster_to_web.py --input ../raw/Final_WFR_Masked_25-05-09_11-12.tif \
    --output ../maps/assets/wfr_raster.png --cmap YlOrRd --crs albers
```

To view the map, open `maps/wfr_populated_places_map.html` in a browser.

**Tools:** ArcGIS Pro · Python · NumPy · GeoPandas · Leaflet / leafmap · AHP · Sobol sensitivity analysis
