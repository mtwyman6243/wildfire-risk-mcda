"""
raster_to_web.py
------------------
Takes a NAD83 California Albers (EPSG:3310) GeoTIFF raster (read via
geotiff_reader.py) and produces a web-ready PNG image overlay in WGS84
for use with Leaflet's L.imageOverlay — no GDAL warp required.

Approach: rather than inverse-project every source pixel (slow, iterative),
we build the *output* WGS84 grid and forward-project each output cell back
into Albers space (closed-form, vectorizable) to look up the nearest source
pixel. This is the standard "backward mapping" resampling strategy.

Usage:
    python raster_to_web.py --input path/to/raster.tif \
        --output ../maps/assets/wfr_raster.png --cmap YlOrRd
"""

import argparse
import json
import math

import numpy as np
from PIL import Image

from geotiff_reader import read_geotiff
from gpkg_reader import albers_to_wgs84
from austria_projection import austria_lambert_to_wgs84, lcc_inverse_to_mgi, mgi_to_wgs84, \
    A_BESSEL, E2_BESSEL, E_BESSEL, LAT1 as AT_LAT1, LAT2 as AT_LAT2, LAT0 as AT_LAT0, \
    LON0 as AT_LON0, X0 as AT_X0, Y0 as AT_Y0, _n as AT_N, _F as AT_F, _RHO0 as AT_RHO0
from gpkg_reader import A, E2, E, LAT1, LAT2, LAT0, LON0, X0, Y0, _N, _C, _RHO0


def wgs84_to_albers(lon_deg, lat_deg):
    """Vectorized forward Albers Conic Equal-Area (ellipsoidal), matches
    the inverse used in gpkg_reader.albers_to_wgs84."""
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    sinphi = np.sin(lat)
    q = (1 - E2) * (
        sinphi / (1 - E2 * sinphi ** 2)
        - (1 / (2 * E)) * np.log((1 - E * sinphi) / (1 + E * sinphi))
    )
    theta = _N * (lon - LON0)
    rho = A * np.sqrt(np.maximum(_C - _N * q, 0)) / _N
    x = rho * np.sin(theta) + X0
    y = _RHO0 - rho * np.cos(theta) + Y0
    return x, y


def wgs84_to_austria_lambert(lon_deg, lat_deg):
    """Vectorized forward transform: WGS84 -> MGI (inverse Helmert, approximated
    as a translation-only shift applied in geocentric space, vectorized) ->
    Lambert Conformal Conic (EPSG:31287)."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    sinlat = np.sin(lat)
    N_ = 6378137.0 / np.sqrt(1 - E2 * sinlat ** 2)  # WGS84 e2 reused (E2 imported = Albers's, but WGS84 close enough for h=0 forward geocentric)
    # Geocentric on WGS84
    Xg = N_ * np.cos(lat) * np.cos(lon)
    Yg = N_ * np.cos(lat) * np.sin(lon)
    Zg = N_ * (1 - E2) * sinlat
    # reverse Helmert shift (WGS84 -> MGI): subtract the MGI->WGS84 translation
    from austria_projection import DX, DY, DZ
    Xm, Ym, Zm = Xg - DX, Yg - DY, Zg - DZ
    # geocentric (Bessel) -> geodetic (MGI)
    lon_m = np.arctan2(Ym, Xm)
    p = np.sqrt(Xm ** 2 + Ym ** 2)
    lat_m = np.arctan2(Zm, p * (1 - E2_BESSEL))
    for _ in range(6):
        sinlat_m = np.sin(lat_m)
        Nb = A_BESSEL / np.sqrt(1 - E2_BESSEL * sinlat_m ** 2)
        lat_m = np.arctan2(Zm + E2_BESSEL * Nb * sinlat_m, p)
    # forward LCC (ellipsoidal, Bessel)
    t = np.tan(np.pi / 4 - lat_m / 2) / (
        ((1 - E_BESSEL * np.sin(lat_m)) / (1 + E_BESSEL * np.sin(lat_m))) ** (E_BESSEL / 2)
    )
    rho = A_BESSEL * AT_F * t ** AT_N
    theta = AT_N * (lon_m - AT_LON0)
    x = rho * np.sin(theta) + AT_X0
    y = AT_RHO0 - rho * np.cos(theta) + AT_Y0
    return x, y


PROJECTIONS = {
    "albers": {"inverse": albers_to_wgs84, "forward": wgs84_to_albers},
    "austria": {"inverse": austria_lambert_to_wgs84, "forward": wgs84_to_austria_lambert},
}


def raster_wgs84_bounds(meta, inverse_fn):
    """Reproject the 4 corners (+ edge midpoints for a tighter box) to find
    a WGS84 bounding box covering the raster."""
    w, h = meta["width"], meta["height"]
    ox, oy, pw, ph = meta["origin_x"], meta["origin_y"], meta["pixel_w"], meta["pixel_h"]
    corners = []
    for fx in (0, 0.5, 1):
        for fy in (0, 0.5, 1):
            x = ox + fx * w * pw
            y = oy - fy * h * ph
            corners.append(inverse_fn(x, y))
    lons = [c[0] for c in corners]
    lats = [c[1] for c in corners]
    return min(lons), min(lats), max(lons), max(lats)


def render(input_path, output_path, cmap_name="YlOrRd", out_width=1400, vmin=None, vmax=None, crs="albers"):
    proj = PROJECTIONS[crs]
    inverse_fn, forward_fn = proj["inverse"], proj["forward"]

    meta = read_geotiff(input_path)
    arr = meta["array"]
    nodata = meta["nodata"]
    mask_valid = np.isfinite(arr) & (arr != nodata)

    if vmin is None:
        vmin = arr[mask_valid].min()
    if vmax is None:
        vmax = arr[mask_valid].max()

    min_lon, min_lat, max_lon, max_lat = raster_wgs84_bounds(meta, inverse_fn)
    aspect = (max_lat - min_lat) / (max_lon - min_lon)
    out_height = max(1, int(out_width * aspect))

    lon_grid, lat_grid = np.meshgrid(
        np.linspace(min_lon, max_lon, out_width),
        np.linspace(max_lat, min_lat, out_height),  # top row = max_lat
    )
    x_alb, y_alb = forward_fn(lon_grid, lat_grid)

    col = ((x_alb - meta["origin_x"]) / meta["pixel_w"]).astype(int)
    row = ((meta["origin_y"] - y_alb) / meta["pixel_h"]).astype(int)

    in_bounds = (col >= 0) & (col < meta["width"]) & (row >= 0) & (row < meta["height"])
    col_c = np.clip(col, 0, meta["width"] - 1)
    row_c = np.clip(row, 0, meta["height"] - 1)

    sampled = arr[row_c, col_c]
    valid = in_bounds & np.isfinite(sampled) & (sampled != nodata)

    from matplotlib import cm
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    colormap = cm.get_cmap(cmap_name)
    rgba = colormap(norm(sampled))
    rgba = (rgba * 255).astype(np.uint8)
    rgba[..., 3] = np.where(valid, 190, 0)  # alpha

    img = Image.fromarray(rgba, mode="RGBA")
    img.save(output_path)

    bounds = {"south": min_lat, "west": min_lon, "north": max_lat, "east": max_lon,
              "vmin": float(vmin), "vmax": float(vmax)}
    return bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--cmap", default="YlOrRd")
    ap.add_argument("--width", type=int, default=1400)
    ap.add_argument("--crs", default="albers", choices=["albers", "austria"])
    args = ap.parse_args()

    bounds = render(args.input, args.output, cmap_name=args.cmap, out_width=args.width, crs=args.crs)
    bounds_path = args.output.rsplit(".", 1)[0] + "_bounds.json"
    with open(bounds_path, "w") as f:
        json.dump(bounds, f)
    print(f"Wrote {args.output}")
    print(f"Wrote {bounds_path}: {bounds}")


if __name__ == "__main__":
    main()
