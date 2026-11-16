"""
Minimal pure-Python GeoPackage (.gpkg) reader — no GDAL/geopandas required.
Reads GPKG binary geometry blobs (WKB with GPKG header) into plain
lon/lat coordinate lists, reprojecting from NAD83 California Albers
(EPSG:3310) to WGS84 using a hand-rolled ellipsoidal Albers inverse
projection (Snyder 1987 formulas) since pyproj isn't available offline.
"""
import sqlite3
import struct
import math

# ---------------------------------------------------------------------
# EPSG:3310 NAD83 / California Albers -> WGS84 inverse projection
# ---------------------------------------------------------------------
A = 6378137.0
F = 1 / 298.257222101  # GRS80
E2 = 2 * F - F * F
E = math.sqrt(E2)
LAT1 = math.radians(34.0)
LAT2 = math.radians(40.5)
LAT0 = math.radians(0.0)
LON0 = math.radians(-120.0)
X0 = 0.0
Y0 = -4000000.0


def _q(phi):
    sinphi = math.sin(phi)
    return (1 - E2) * (
        sinphi / (1 - E2 * sinphi ** 2)
        - (1 / (2 * E)) * math.log((1 - E * sinphi) / (1 + E * sinphi))
    )


def _m(phi):
    sinphi = math.sin(phi)
    return math.cos(phi) / math.sqrt(1 - E2 * sinphi ** 2)


_m1, _m2 = _m(LAT1), _m(LAT2)
_q0, _q1, _q2 = _q(LAT0), _q(LAT1), _q(LAT2)
_N = (_m1 ** 2 - _m2 ** 2) / (_q2 - _q1)
_C = _m1 ** 2 + _N * _q1
_RHO0 = A * math.sqrt(_C - _N * _q0) / _N


def albers_to_wgs84(x, y):
    xp = x - X0
    yp = _RHO0 - (y - Y0)
    rho = math.sqrt(xp ** 2 + yp ** 2)
    if _N < 0:
        rho = -rho
        xp, yp = -xp, -yp
    theta = math.atan2(xp, yp)
    q = (_C - (rho * _N / A) ** 2) / _N
    q = max(-2.0, min(2.0, q))  # valid domain for asin(q/2)
    phi = math.asin(q / 2)
    for _ in range(6):
        sinphi = math.sin(phi)
        denom = 1 - E2 * sinphi ** 2
        phi += (denom ** 2 / (2 * math.cos(phi))) * (
            q / (1 - E2)
            - sinphi / denom
            + (1 / (2 * E)) * math.log((1 - E * sinphi) / (1 + E * sinphi))
        )
    lon = LON0 + theta / _N
    return math.degrees(lon), math.degrees(phi)


# ---------------------------------------------------------------------
# GPKG binary geometry parsing
# ---------------------------------------------------------------------
_ENV_SIZES = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}


def _read_wkb_geom(buf, offset, byteorder):
    fmt_prefix = '<' if byteorder == 1 else '>'
    gtype = struct.unpack(fmt_prefix + 'I', buf[offset:offset + 4])[0]
    offset += 4
    base_type = gtype % 1000  # strip Z/M/ZM flags (1000/2000/3000)
    has_z = gtype >= 1000

    def read_point():
        nonlocal offset
        n = 3 if has_z else 2
        vals = struct.unpack(fmt_prefix + 'd' * n, buf[offset:offset + 8 * n])
        offset += 8 * n
        return vals[0], vals[1]

    def read_ring():
        nonlocal offset
        npts = struct.unpack(fmt_prefix + 'I', buf[offset:offset + 4])[0]
        offset += 4
        pts = [read_point() for _ in range(npts)]
        return pts

    def read_polygon():
        nonlocal offset
        nrings = struct.unpack(fmt_prefix + 'I', buf[offset:offset + 4])[0]
        offset += 4
        return [read_ring() for _ in range(nrings)]

    if base_type == 3:  # Polygon
        poly = read_polygon()
        return [poly], offset
    elif base_type == 6:  # MultiPolygon
        npoly = struct.unpack(fmt_prefix + 'I', buf[offset:offset + 4])[0]
        offset += 4
        polys = []
        for _ in range(npoly):
            sub_order = buf[offset]
            sub_polys, offset_after = _read_wkb_geom(buf, offset + 1, sub_order)
            polys.extend(sub_polys)
            offset = offset_after
        return polys, offset
    else:
        raise NotImplementedError(f"geometry base type {base_type} not handled")


def parse_gpkg_geometry(blob, project_fn=albers_to_wgs84):
    """Returns list of polygons; each polygon is a list of rings; each ring
    is a list of (lon, lat) tuples, reprojected to WGS84 via project_fn."""
    flags = blob[3]
    envelope_code = (flags >> 1) & 0b111
    env_len = _ENV_SIZES[envelope_code]
    wkb_start = 8 + env_len
    byteorder = blob[wkb_start]
    polys, _ = _read_wkb_geom(buf=blob, offset=wkb_start + 1, byteorder=byteorder)
    out = []
    for poly in polys:
        out_poly = []
        for ring in poly:
            out_poly.append([project_fn(x, y) for x, y in ring])
        out.append(out_poly)
    return out


def read_gpkg_layer(path, project_fn=albers_to_wgs84):
    """Reads the single feature layer in a GeoPackage into
    (attributes_list, geometry_list) where geometry_list holds
    lon/lat MultiPolygon rings per feature."""
    con = sqlite3.connect(path)
    cur = con.cursor()
    table = cur.execute(
        "SELECT table_name FROM gpkg_contents WHERE data_type='features' LIMIT 1"
    ).fetchone()[0]
    cols = [r[1] for r in cur.execute(f'PRAGMA table_info("{table}")')]
    attr_cols = [c for c in cols if c != 'geom']
    rows = cur.execute(
        f'SELECT {", ".join(f"""\"{c}\"""" for c in attr_cols)}, geom FROM "{table}"'
    ).fetchall()
    features = []
    for row in rows:
        attrs = dict(zip(attr_cols, row[:-1]))
        geom_blob = row[-1]
        multipoly = parse_gpkg_geometry(geom_blob, project_fn=project_fn)
        features.append((attrs, multipoly))
    con.close()
    return features
