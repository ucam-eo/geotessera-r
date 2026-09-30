"""Build the synthetic example store and golden outputs from the Python reader.

The store ships as inst/extdata/example.zarr for examples and tests.

Run from a geotessera checkout, which provides zarr, pyproj and geotessera:

    uv run --project ../geotessera python data-raw/make-fixture.py
"""

import json
import math
import shutil
from pathlib import Path

import numpy as np
import zarr
from pyproj import Transformer
from zarr.codecs import BloscCodec

from geotessera.store import GeoTesseraZarr

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "tests" / "testthat" / "fixtures"
STORE = ROOT / "inst" / "extdata" / "example.zarr"
YEARS = [2023, 2024]
BANDS = 4
DEPTHS = {2: "embeddings_d2", 4: "embeddings"}
PX = 10.0
H = W = 256
SHARD = 128
LAT = 52.0

EMB_BLOSC = BloscCodec(cname="zstd", clevel=3, shuffle="bitshuffle", typesize=1)
FLOAT_BLOSC = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle", typesize=4)
COORD_BLOSC = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle", typesize=8)


def zone_grid(zone):
    """Zone 30 ends 8 px east of the lon 0 seam; zone 31 starts on it."""
    epsg = 32600 + zone
    e, n = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True).transform(
        0.0, LAT
    )
    if zone == 30:
        ox = math.ceil(e / PX) * PX + 8 * PX - W * PX
    else:
        ox = math.floor(e / PX) * PX
    oy = round(n / PX) * PX + (H // 2) * PX
    return epsg, ox, oy


def zone_values(zone):
    t, b, y, x = np.meshgrid(
        np.arange(len(YEARS)), np.arange(BANDS), np.arange(H), np.arange(W), indexing="ij"
    )
    emb = ((x * 7 + y * 3 + b * 11 + t * 5 + zone) % 251 - 125).astype(np.int8)
    yy, xx = np.meshgrid(np.arange(H), np.arange(W), indexing="ij")
    base = 0.05 if zone == 30 else 0.07
    scales = np.empty((len(YEARS), H, W), np.float32)
    for ti in range(len(YEARS)):
        scales[ti] = base * (1 + ((xx + yy + ti) % 5) / 10)
    scales[:, 40:60, 40:60] = np.nan  # water
    scales[:, 100, 100] = np.inf  # isolated hole
    scales[:, 150:153, 150:153] = np.inf  # hole wider than the search radius
    if zone == 31:
        scales[1, 192:208, 16:32] = np.inf  # absent inner chunk in 2024
        emb[1, :, 192:208, 16:32] = 0
        scales[1, :, 0:4] = np.inf  # unwritten west edge in 2024
        emb[1, :, :, 0:4] = 0
    return emb, scales


def sharded(group, name, data, chunks, shards, compressor, fill, dims):
    arr = group.create_array(
        name,
        shape=data.shape,
        dtype=data.dtype,
        chunks=chunks,
        shards=shards,
        compressors=compressor,
        fill_value=fill,
        dimension_names=dims,
    )
    arr[:] = data


def build_store():
    shutil.rmtree(STORE, ignore_errors=True)
    root = zarr.open_group(STORE, mode="w", zarr_format=3)
    root.attrs.update(
        {
            "geoemb:type": "pixel",
            "geoemb:dimensions": BANDS,
            "geoemb:model": "https://geotessera.org/model/test",
            "geoemb:build_version": "test",
            "geoemb:depths": [{"dimensions": d, "array": a} for d, a in DEPTHS.items()],
        }
    )
    for zone in (30, 31):
        epsg, ox, oy = zone_grid(zone)
        g = root.create_group(f"utm{zone:02d}")
        g.attrs.update(
            {
                "spatial:dimensions": ["y", "x"],
                "spatial:transform_type": "affine",
                "spatial:transform": [PX, 0.0, ox, 0.0, -PX, oy],
                "spatial:shape": [H, W],
                "spatial:registration": "pixel",
                "proj:code": f"EPSG:{epsg}",
            }
        )
        emb, scales = zone_values(zone)
        for depth, name in DEPTHS.items():
            inner = 16 if depth == BANDS else 32
            sharded(
                g, name, emb[:, :depth], (1, depth, inner, inner),
                (1, depth, SHARD, SHARD), EMB_BLOSC, 0,
                ("time", "band" if depth == BANDS else f"band_d{depth}", "y", "x"),
            )
        sharded(g, "scales", scales, (1, 16, 16), (1, SHARD, SHARD), FLOAT_BLOSC,
                np.float32("inf"), ("time", "y", "x"))
        for name, values, dtype in (
            ("time", YEARS, np.int32),
            ("band", range(BANDS), np.int32),
            ("x", ox + (np.arange(W) + 0.5) * PX, np.float64),
            ("y", oy - (np.arange(H) + 0.5) * PX, np.float64),
        ):
            values = np.asarray(values, dtype=dtype)
            a = g.create_array(name, shape=values.shape, dtype=dtype,
                               compressors=COORD_BLOSC if dtype == np.float64 else FLOAT_BLOSC,
                               dimension_names=(name,))
            a[:] = values
    zarr.consolidate_metadata(STORE)


def lonlat(zone, col, row):
    epsg, ox, oy = zone_grid(zone)
    to_geo = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
    return to_geo.transform(ox + (col + 0.5) * PX, oy - (row + 0.5) * PX)


def jsonable(a):
    a = np.asarray(a, dtype=float)
    return np.where(np.isnan(a), None, a).tolist()


def golden():
    gt = GeoTesseraZarr(str(STORE))
    points = {
        "valid": lonlat(31, 60, 80),
        "water": lonlat(31, 50, 50),
        "hole": lonlat(31, 100, 100),
        "wide_hole": lonlat(31, 151, 151),
        "absent_chunk": lonlat(31, 20, 200),
        "zone30": lonlat(30, 120, 30),
        "seam_edge": lonlat(31, 1, 128),
        "outside": (0.5, 52.0),
    }
    names = list(points)
    coords = [points[k] for k in names]
    out = {"points": dict(zip(names, map(list, coords))), "years": gt.years}
    for year in YEARS:
        out[f"sample_{year}"] = jsonable(gt.sample_points(coords, year))
        out[f"probe_{year}"] = [gt.probe(lon, lat, year)[1] for lon, lat in coords]
        out[f"sample_nosearch_{year}"] = jsonable(
            gt.sample_points(coords, year, search_px=0, cross_zone=False)
        )
    out["sample_d2"] = jsonable(gt.sample_points(coords, 2024, depth=2))

    # Bounds between pixel centres, so projection round-off cannot move them.
    lon0, lat0 = lonlat(31, 29.75, 29.75)
    lon1, lat1 = lonlat(31, 90.25, 70.25)
    bbox = (lon0, lat1, lon1, lat0)
    mosaic, transform, crs = gt.read_region(bbox, 2024)
    emb, scales, _, _ = gt.read_region_quantized(bbox, 2024)
    strips = [(b.shape, list(t)[:6]) for b, t, _ in gt.iter_region(bbox, 2024, strip_rows=16)]
    out["region"] = {
        "bbox": bbox,
        "shape": list(mosaic.shape),
        "transform": list(transform)[:6],
        "crs": crs,
        "values": jsonable(mosaic),
        "quantized": np.asarray(emb, dtype=int).tolist(),
        "scales": jsonable(scales),
        "strips": [[list(s), t] for s, t in strips],
    }

    for key, (lon, lat, size) in {
        "patch_native": (*lonlat(31, 120, 120), 16),
        "patch_edge": (*lonlat(31, 250, 5), 16),
        "patch_seam": (0.0, LAT, 24),
    }.items():
        patch, transform, crs = gt.read_patch(lon, lat, 2023, size)
        out[key] = {
            "lon": lon, "lat": lat, "size": size,
            "transform": list(transform)[:6],
            "native": crs.startswith("EPSG:"),
            "crs": crs,
            "values": jsonable(patch),
        }

    west, _ = lonlat(30, 200.25, 100)
    east, _ = lonlat(31, 60.25, 100)
    _, north = lonlat(31, 60, 90.25)
    _, south = lonlat(31, 60, 110.25)
    ebbox = (west, south, east, north)
    out["export"] = {
        "bbox": ebbox,
        "estimates": gt.export_geotiffs(ebbox, 2024, OUT / "unused", dry_run=True),
    }
    (OUT / "golden.json").write_text(json.dumps(out))


if __name__ == "__main__":
    build_store()
    golden()
