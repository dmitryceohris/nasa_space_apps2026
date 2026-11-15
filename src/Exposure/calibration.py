"""Explicit, versioned detector calibration; never infer a calibration from a date."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import json
import os
from functools import lru_cache

import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates


def parent_pixel(exp, x, y):
    """IRSA alternate A coordinates map cutout pixels onto the parent detector (0 based)."""
    h = exp.header
    if h is not None and "CRPIX1A" in h and "CRPIX2A" in h:
        return np.asarray(x) + 1 - h["CRPIX1A"], np.asarray(y) + 1 - h["CRPIX2A"]
    if exp.shape == (2040, 2040):
        return np.asarray(x), np.asarray(y)
    raise ValueError("parent_detector_mapping_missing")


def sample(array, x, y):
    x, y = np.broadcast_arrays(x, y)
    values = map_coordinates(np.asarray(array, float), np.array([y.ravel(), x.ravel()]),
                             order=1, mode="constant", cval=np.nan)
    return values.reshape(x.shape)


@dataclass
class Calibration:
    cwave: np.ndarray
    cband: np.ndarray
    detector: int
    solid_angle: np.ndarray | None = None  # arcsec^2
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.cwave.shape != self.cband.shape or self.cwave.ndim != 2:
            raise ValueError("CWAVE/CBAND shape mismatch")
        if self.solid_angle is not None and self.solid_angle.shape != self.cwave.shape:
            raise ValueError("solid-angle shape mismatch")

    def wavelength_at(self, exp, x, y):
        if exp.detector != self.detector:
            raise ValueError("calibration_detector_mismatch")
        px, py = parent_pixel(exp, x, y)
        return float(sample(self.cwave, px, py)), float(sample(self.cband, px, py))

    def omega_sr(self, exp, x, y):
        if self.solid_angle is None:
            raise ValueError("missing_solid_angle")
        from astropy import units as u
        px, py = parent_pixel(exp, x, y)
        return sample(self.solid_angle, px, py) * u.arcsec.to(u.rad) ** 2


def load_calibration(wave_path, solid_angle_path=None, *, detector: int, expected_sha256=None):
    paths = [Path(wave_path)] + ([Path(solid_angle_path)] if solid_angle_path else [])
    records = []
    for p in paths:
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        if expected_sha256 and expected_sha256.get(p.name) != digest:
            raise ValueError(f"calibration_checksum_mismatch: {p.name}")
        records.append({"file": p.name, "sha256": digest})
    with fits.open(paths[0], memmap=False) as h:
        for hd in h:
            if "DETECTOR" in hd.header and int(hd.header["DETECTOR"]) != detector:
                raise ValueError("calibration_detector_mismatch")
        wave, band = h["CWAVE"].data.copy(), h["CBAND"].data.copy()
    angle = None
    if solid_angle_path:
        with fits.open(solid_angle_path, memmap=False) as h:
            hd = next(hd for hd in h if hd.data is not None and hd.data.ndim == 2)
            if "DETECTOR" in hd.header and int(hd.header["DETECTOR"]) != detector:
                raise ValueError("calibration_detector_mismatch")
            unit = hd.header.get("BUNIT", "").lower().replace(" ", "")
            if unit not in ("arcsec2", "arcsec^2", "arcsec**2"):
                raise ValueError("unverified_solid_angle_unit")
            angle = hd.data.copy()
    return Calibration(wave, band, detector, angle, {"products": records})


@lru_cache(maxsize=12)
def _registered_calibration(manifest_path, manifest_mtime, detector, release, pipeline_version):
    path = Path(manifest_path)
    rows = json.loads(path.read_text())["calibrations"]
    candidates = [r for r in rows if r["detector"] == detector and release in r["releases"]
                  and pipeline_version in r["pipeline_versions"]]
    if len(candidates) != 1:
        return None
    row = candidates[0]
    wave, angle = path.parent / row["wave_file"], path.parent / row["solid_angle_file"]
    result = load_calibration(wave, angle, detector=detector,
                              expected_sha256={wave.name:row["wave_sha256"],angle.name:row["solid_angle_sha256"]})
    result.provenance.update({k:row[k] for k in ("wave_url","solid_angle_url","releases","pipeline_versions")})
    return result


def attach_registered_calibration(exp, pipeline_version):
    default = Path(__file__).resolve().parents[2] / "cache/calibrations/manifest.json"
    path = Path(os.environ.get("SPHEREX_CALIBRATION_MANIFEST", default))
    release = exp.provenance.get("release")
    if path.is_file() and release:
        exp.calibration = _registered_calibration(str(path),path.stat().st_mtime_ns,
                                                exp.detector,release,pipeline_version)
    return exp
