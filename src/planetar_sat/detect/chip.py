"""SAR scene → detection batch.

Wraps the CFAR detector with rasterio I/O: opens a Sentinel-1 GRD GeoTIFF,
runs CFAR over the VV pol, and projects pixel coordinates back to lat/lon.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from planetar_sat.detect.cfar import Detection, cfar_ca

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeoDetection:
    scene_id: str
    lat: float
    lon: float
    snr: float
    bbox_px: tuple[int, int, int, int]
    acquired_at_ns: int


def detect_scene(
    tif_path: Path,
    scene_id: str,
    acquired_at_ns: int,
    cfar_kwargs: dict | None = None,
) -> list[GeoDetection]:
    """Run CFAR on a single VV GRD GeoTIFF and return geo-referenced detections."""
    try:
        import rasterio
        from rasterio.transform import xy
    except ImportError as e:  # pragma: no cover
        raise RuntimeError(
            "rasterio is required for scene detection. `pip install rasterio` "
            "or `make install`."
        ) from e

    with rasterio.open(tif_path) as ds:
        band = ds.read(1).astype(np.float32)
        transform = ds.transform
        crs = ds.crs

    dets = cfar_ca(band, **(cfar_kwargs or {}))
    geo_dets: list[GeoDetection] = []
    for d in dets:
        # rasterio xy() returns (x, y) in dataset CRS. For lat/lon we'd need a
        # reproject if CRS isn't already 4326. Sentinel-1 GRD on Copernicus is
        # typically WGS84-aligned EPSG:4326 after Level-1 GRD processing, but
        # full IW SLC products use UTM. We project explicitly when needed.
        x, y = xy(transform, d.row, d.col)
        if crs is not None and crs.to_epsg() != 4326:
            try:
                from rasterio.warp import transform as warp_transform
                xs, ys = warp_transform(crs, "EPSG:4326", [x], [y])
                lon, lat = xs[0], ys[0]
            except Exception:  # pragma: no cover
                lon, lat = float("nan"), float("nan")
        else:
            lon, lat = x, y
        geo_dets.append(
            GeoDetection(
                scene_id=scene_id,
                lat=float(lat),
                lon=float(lon),
                snr=d.snr,
                bbox_px=d.bbox,
                acquired_at_ns=acquired_at_ns,
            )
        )
    log.info("scene %s: %d geo-detections", scene_id, len(geo_dets))
    return geo_dets
