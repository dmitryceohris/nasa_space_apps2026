"""Explicit, versioned detector calibration; never infer a calibration from a date."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
import json
import os
from functools import lru_cache
import threading
import io

from runtime_storage import current_workspace

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


def _wave_arrays(source, detector):
    with fits.open(source, memmap=False) as h:
        for hd in h:
            if "DETECTOR" in hd.header and int(hd.header["DETECTOR"]) != detector:
                raise ValueError("calibration_detector_mismatch")
        return h["CWAVE"].data.copy(), h["CBAND"].data.copy()


def _solid_angle_array(source, detector):
    with fits.open(source, memmap=False) as h:
        hd = next(hd for hd in h if hd.data is not None and hd.data.ndim == 2)
        if "DETECTOR" in hd.header and int(hd.header["DETECTOR"]) != detector:
            raise ValueError("calibration_detector_mismatch")
        unit = hd.header.get("BUNIT", "").lower().replace(" ", "")
        if unit not in ("arcsec2", "arcsec^2", "arcsec**2"):
            raise ValueError("unverified_solid_angle_unit")
        return hd.data.copy()


def load_calibration(wave_path, solid_angle_path=None, *, detector: int, expected_sha256=None):
    paths = [Path(wave_path)] + ([Path(solid_angle_path)] if solid_angle_path else [])
    records = []
    for p in paths:
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        if expected_sha256 and expected_sha256.get(p.name) != digest:
            raise ValueError(f"calibration_checksum_mismatch: {p.name}")
        records.append({"file": p.name, "sha256": digest})
    wave, band = _wave_arrays(paths[0], detector)
    angle = _solid_angle_array(paths[1], detector) if solid_angle_path else None
    return Calibration(wave, band, detector, angle, {"products": records})


_calibration_lock = threading.RLock()


def _registered_calibration(manifest_path, manifest_mtime, detector, release, pipeline_version):
    # functools.lru_cache alone allows duplicate concurrent misses and downloads.
    with _calibration_lock:
        return _load_registered_calibration(manifest_path, manifest_mtime, detector, release, pipeline_version)


@lru_cache(maxsize=12)
def _load_registered_calibration(manifest_path, manifest_mtime, detector, release, pipeline_version):
    path = Path(manifest_path)
    rows = json.loads(path.read_text())["calibrations"]
    candidates = [r for r in rows if r["detector"] == detector and release in r["releases"]
                  and pipeline_version in r["pipeline_versions"]]
    if len(candidates) != 1:
        return None
    row = candidates[0]
    wave, angle = path.parent / row["wave_file"], path.parent / row["solid_angle_file"]
    if os.environ.get("SPHEREX_CALIBRATION_CACHE"):
        import tempfile
        import httpx
        directory = Path(os.environ["SPHEREX_CALIBRATION_CACHE"])
        directory.mkdir(exist_ok=True, parents=True)
        wave, angle = directory / wave.name, directory / angle.name
    workspace = current_workspace()
    if workspace is not None:
        # Vercel's hydrated dependencies already occupy most of /tmp. Keep the verified
        # calibration in the existing in-memory LRU rather than duplicating it on disk.
        wave, band = _wave_arrays(io.BytesIO(workspace.verified_bytes(row["wave_url"], row["wave_sha256"])), detector)
        angle = _solid_angle_array(io.BytesIO(workspace.verified_bytes(row["solid_angle_url"], row["solid_angle_sha256"])), detector)
        records = [{"file":row[kind+"_file"],"sha256":row[kind+"_sha256"]} for kind in ("wave","solid_angle")]
        result = Calibration(wave, band, detector, angle, {"products":records})
        result.provenance.update({k:row[k] for k in ("wave_url","solid_angle_url","releases","pipeline_versions")})
        return result
    if os.environ.get("SPHEREX_CALIBRATION_CACHE"):
        for kind, dest in (("wave", wave), ("solid_angle", angle)):
            expected = row[f"{kind}_sha256"]
            if dest.is_file() and hashlib.sha256(dest.read_bytes()).hexdigest() == expected:
                os.utime(dest, None)
                continue
            dest.unlink(missing_ok=True)
            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as stream:
                temporary = Path(stream.name)
            try:
                with httpx.stream("GET", row[f"{kind}_url"], timeout=45) as response:
                    response.raise_for_status()
                    with temporary.open("wb") as output:
                        for chunk in response.iter_bytes():
                            output.write(chunk)
                if hashlib.sha256(temporary.read_bytes()).hexdigest() != expected:
                    raise ValueError("calibration_checksum_mismatch")
                temporary.replace(dest)
            finally:
                temporary.unlink(missing_ok=True)
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
