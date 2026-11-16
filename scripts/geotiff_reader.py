"""
geotiff_reader.py
------------------
Minimal GeoTIFF reader for single-band, uncompressed rasters (no
GDAL/rasterio required) — reads the IFD tags with `struct`, then pulls
pixel data with `numpy` for speed. Handles the strip-based layout these
WFR rasters use (RowsPerStrip=1, SampleFormat=IEEE float, 64 bits).

This covers the common case needed here. For general-purpose GeoTIFF
reading (compressed, tiled, multi-band), use rioxarray/rasterio instead
— see build_project_map_template.py for that pattern.
"""

import struct
import numpy as np

TYPE_SIZES = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8}

TAG_IMAGEWIDTH = 256
TAG_IMAGELENGTH = 257
TAG_BITSPERSAMPLE = 258
TAG_COMPRESSION = 259
TAG_STRIPOFFSETS = 273
TAG_ROWSPERSTRIP = 278
TAG_STRIPBYTECOUNTS = 279
TAG_SAMPLEFORMAT = 339
TAG_MODELPIXELSCALE = 33550
TAG_MODELTIEPOINT = 33922
TAG_GDAL_NODATA = 42113


def _read_ifd(f, byteorder):
    num_entries = struct.unpack(byteorder + 'H', f.read(2))[0]
    tags = {}
    for _ in range(num_entries):
        entry = f.read(12)
        tag, typ, count = struct.unpack(byteorder + 'HHI', entry[0:8])
        valbytes = entry[8:12]
        tsize = TYPE_SIZES.get(typ, 1)
        total = tsize * count
        if total <= 4:
            data = valbytes[:total]
        else:
            offset = struct.unpack(byteorder + 'I', valbytes)[0]
            cur = f.tell()
            f.seek(offset)
            data = f.read(total)
            f.seek(cur)
        tags[tag] = (typ, count, data)
    return tags


def _unpack_array(byteorder, typ, count, data):
    fmt = {1: 'B', 3: 'H', 4: 'I', 11: 'f', 12: 'd'}[typ]
    return struct.unpack(byteorder + fmt * count, data)


def read_geotiff(path):
    """Returns dict with keys: array (2D numpy float64), width, height,
    origin_x, origin_y, pixel_w, pixel_h, nodata."""
    with open(path, 'rb') as f:
        header = f.read(8)
        byteorder = '<' if header[0:2] == b'II' else '>'
        ifd_offset = struct.unpack(byteorder + 'I', header[4:8])[0]
        f.seek(ifd_offset)
        tags = _read_ifd(f, byteorder)

        width = _unpack_array(byteorder, 3, 1, tags[TAG_IMAGEWIDTH][2])[0]
        height = _unpack_array(byteorder, 3, 1, tags[TAG_IMAGELENGTH][2])[0]
        bits = _unpack_array(byteorder, 3, 1, tags[TAG_BITSPERSAMPLE][2])[0]
        sample_format = _unpack_array(byteorder, 3, 1, tags[TAG_SAMPLEFORMAT][2])[0] if TAG_SAMPLEFORMAT in tags else 1
        rows_per_strip = _unpack_array(byteorder, 3, 1, tags[TAG_ROWSPERSTRIP][2])[0]

        so_typ, so_count, so_data = tags[TAG_STRIPOFFSETS]
        sc_typ, sc_count, sc_data = tags[TAG_STRIPBYTECOUNTS]
        strip_offsets = _unpack_array(byteorder, so_typ, so_count, so_data)
        strip_bytecounts = _unpack_array(byteorder, sc_typ, sc_count, sc_data)

        dtype_map = {(64, 3): np.float64, (32, 3): np.float32,
                     (16, 1): np.uint16, (8, 1): np.uint8, (32, 1): np.uint32}
        dtype = dtype_map.get((bits, sample_format), np.float64)
        np_byteorder = '<' if byteorder == '<' else '>'
        dtype = np.dtype(dtype).newbyteorder(np_byteorder)

        rows = []
        for off, cnt in zip(strip_offsets, strip_bytecounts):
            f.seek(off)
            buf = f.read(cnt)
            rows.append(np.frombuffer(buf, dtype=dtype))
        array = np.concatenate(rows).reshape(height, width).astype(np.float64)

        scale = _unpack_array(byteorder, 12, 3, tags[TAG_MODELPIXELSCALE][2])
        tiepoint = _unpack_array(byteorder, 12, 6, tags[TAG_MODELTIEPOINT][2])
        origin_x, origin_y = tiepoint[3], tiepoint[4]
        pixel_w, pixel_h = scale[0], scale[1]

        nodata = None
        if TAG_GDAL_NODATA in tags:
            raw = tags[TAG_GDAL_NODATA][2]
            try:
                nodata = float(raw.split(b'\x00')[0])
            except ValueError:
                nodata = None

        return {
            "array": array, "width": width, "height": height,
            "origin_x": origin_x, "origin_y": origin_y,
            "pixel_w": pixel_w, "pixel_h": pixel_h, "nodata": nodata,
        }
