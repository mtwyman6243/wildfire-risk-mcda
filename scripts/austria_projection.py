"""
austria_projection.py
-----------------------
Reprojects MGI / Austria Lambert (EPSG:31287) coordinates to
WGS84 lon/lat, in two steps:

  1. Inverse Lambert Conformal Conic (2SP), ellipsoidal, on the Bessel 1841
     ellipsoid (the MGI datum's ellipsoid) -> geographic lat/lon on MGI.
  2. A 3-parameter (translation-only) Helmert datum shift, MGI -> WGS84
     (EPSG transformation 1618, "MGI to WGS 84 (1)"), via geocentric
     Cartesian coordinates.

This is the standard two-step process GDAL/PROJ do internally for this
EPSG code; done by hand here since pyproj isn't available offline.
"""

import math

# ---------------------------------------------------------------------
# EPSG:31287 (MGI / Austria Lambert) projection parameters
# ---------------------------------------------------------------------
A_BESSEL = 6377397.155
F_BESSEL = 1 / 299.152812800003
E2_BESSEL = 2 * F_BESSEL - F_BESSEL ** 2
E_BESSEL = math.sqrt(E2_BESSEL)

LAT1 = math.radians(49.0)
LAT2 = math.radians(46.0)
LAT0 = math.radians(47.5)
LON0 = math.radians(13.3333333333333)
X0 = 400000.0
Y0 = 400000.0


def _m(phi, e):
    return math.cos(phi) / math.sqrt(1 - e ** 2 * math.sin(phi) ** 2)


def _t(phi, e):
    sinphi = math.sin(phi)
    return math.tan(math.pi / 4 - phi / 2) / (
        ((1 - e * sinphi) / (1 + e * sinphi)) ** (e / 2)
    )


_m1 = _m(LAT1, E_BESSEL)
_m2 = _m(LAT2, E_BESSEL)
_t0 = _t(LAT0, E_BESSEL)
_t1 = _t(LAT1, E_BESSEL)
_t2 = _t(LAT2, E_BESSEL)
_n = (math.log(_m1) - math.log(_m2)) / (math.log(_t1) - math.log(_t2))
_F = _m1 / (_n * _t1 ** _n)
_RHO0 = A_BESSEL * _F * _t0 ** _n


def lcc_inverse_to_mgi(x, y):
    """EPSG:31287 (x,y) -> geographic lat/lon (degrees) on the MGI datum."""
    xp = x - X0
    yp = _RHO0 - (y - Y0)
    rho = math.copysign(math.sqrt(xp ** 2 + yp ** 2), _n)
    theta = math.atan2(xp, yp) if _n > 0 else math.atan2(-xp, -yp)
    t = (rho / (A_BESSEL * _F)) ** (1 / _n)
    phi = math.pi / 2 - 2 * math.atan(t)
    for _ in range(6):
        sinphi = math.sin(phi)
        phi = math.pi / 2 - 2 * math.atan(
            t * (((1 - E_BESSEL * sinphi) / (1 + E_BESSEL * sinphi)) ** (E_BESSEL / 2))
        )
    lon = theta / _n + LON0
    return math.degrees(lon), math.degrees(phi)


# ---------------------------------------------------------------------
# MGI -> WGS84 datum shift (EPSG transformation 1618, 3-parameter)
# ---------------------------------------------------------------------
DX, DY, DZ = 577.326, 90.129, 463.919  # metres, MGI -> WGS84

A_WGS84 = 6378137.0
F_WGS84 = 1 / 298.257223563
E2_WGS84 = 2 * F_WGS84 - F_WGS84 ** 2


def _geodetic_to_geocentric(lat_deg, lon_deg, h, a, e2):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sinlat = math.sin(lat)
    N = a / math.sqrt(1 - e2 * sinlat ** 2)
    X = (N + h) * math.cos(lat) * math.cos(lon)
    Y = (N + h) * math.cos(lat) * math.sin(lon)
    Z = (N * (1 - e2) + h) * sinlat
    return X, Y, Z


def _geocentric_to_geodetic(X, Y, Z, a, e2):
    lon = math.atan2(Y, X)
    p = math.sqrt(X ** 2 + Y ** 2)
    lat = math.atan2(Z, p * (1 - e2))
    for _ in range(6):
        sinlat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sinlat ** 2)
        lat = math.atan2(Z + e2 * N * sinlat, p)
    return math.degrees(lat), math.degrees(lon)


def mgi_to_wgs84(lon_deg, lat_deg):
    X, Y, Z = _geodetic_to_geocentric(lat_deg, lon_deg, 0.0, A_BESSEL, E2_BESSEL)
    X, Y, Z = X + DX, Y + DY, Z + DZ
    lat_w, lon_w = _geocentric_to_geodetic(X, Y, Z, A_WGS84, E2_WGS84)
    return lon_w, lat_w


def austria_lambert_to_wgs84(x, y):
    """EPSG:31287 (x, y) -> WGS84 (lon, lat) in degrees."""
    lon_mgi, lat_mgi = lcc_inverse_to_mgi(x, y)
    return mgi_to_wgs84(lon_mgi, lat_mgi)
