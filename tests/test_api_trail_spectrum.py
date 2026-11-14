"""Offline end-to-end tests for trail mode, /api/spectrum and /api/lightcurve.

IRSA is replaced by synthetic SPHEREx-like FITS files (IMAGE/FLAGS/VARIANCE/ZODI + header keys), so the real
code paths run: job queue -> cutout -> reprojection/photometry (src/processing) -> PNG/JSON.
"""
import os
import time

os.environ["PREWARM"] = "0"

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from fastapi.testclient import TestClient
from PIL import Image

from api import cache, irsa
from api.main import app
from Exposure import DETECTOR_BANDS

RA, DEC = 150.0, 20.0
EPOCH_MJD = [60900.0, 61080.0, 61260.0]   # three ~6-month epochs
STEP_ARCSEC = 18.0                         # star moves ~3 px north per epoch
c = TestClient(app)


def _fits(path, det, epoch, k, flux=60.0, n=41):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA, DEC]
    w.wcs.crpix = [(n + 1) / 2 + 0.3 * k, (n + 1) / 2 - 0.2 * k]   # slightly different pointing per exposure
    w.wcs.cdelt = [-6.15 / 3600, 6.15 / 3600]
    hdr = w.to_header()
    hdr["DETECTOR"], hdr["OBSID"] = det, f"2025W{10 + epoch}_1A_000{k}_1"
    hdr["MJD-AVG"] = EPOCH_MJD[epoch] + k * 0.5
    star = SkyCoord(RA * u.deg, DEC * u.deg).spherical_offsets_by(0 * u.arcsec,
                                                                 (epoch - 1) * STEP_ARCSEC * u.arcsec)
    x0, y0 = w.world_to_pixel(star)
    yy, xx = np.indices((n, n))
    rng = np.random.default_rng(100 * det + 10 * epoch + k)
    img = 1.0 + rng.normal(0, 0.03, (n, n)) + flux * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / 2)
    fh = w.to_header()  # real SPHEREx FLAGS carry the same WCS as IMAGE; the API reprojects FLAGS with it
    fh["MSKN1"], fh["MSKN21"] = "OVERFLOW", "SOURCE"
    fits.HDUList([fits.PrimaryHDU(),
                  fits.ImageHDU(img.astype(np.float32), hdr, name="IMAGE"),
                  fits.ImageHDU(np.zeros((n, n), np.int32), fh, name="FLAGS"),
                  fits.ImageHDU(np.full((n, n), 0.03 ** 2, np.float32), name="VARIANCE"),
                  fits.ImageHDU(np.ones((n, n), np.float32), name="ZODI")]).writeto(path, overwrite=True)


@pytest.fixture()
def fake_irsa(tmp_path, monkeypatch):
    monkeypatch.setattr(cache, "CACHE", tmp_path)
    calls = {"cutout": 0}

    def exposures(ra, dec, det):
        lo, hi, _ = DETECTOR_BANDS[det]
        return [dict(uri=f"syn/{det}/{e}/{k}", mjd=EPOCH_MJD[e] + k * 0.5, lam=[lo, hi])
                for e in range(3) for k in range(2)]

    def cutout(uri, ra, dec, size_deg):
        calls["cutout"] += 1
        _, det, e, k = uri.split("/")
        p = tmp_path / f"syn_{det}_{e}_{k}.fits"
        if not p.exists():
            _fits(p, int(det), int(e), int(k))
        return p

    monkeypatch.setattr(irsa, "exposures", exposures)
    monkeypatch.setattr(irsa, "cutout", cutout)
    return calls


def _wait(url, timeout=60):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = c.get(url)
        if r.status_code != 202:
            return r
        time.sleep(0.2)
    raise AssertionError(f"{url} still running after {timeout}s")


def test_trail_mode(fake_irsa):
    r = _wait(f"/api/view?ra={RA}&dec={DEC}&fov=4&det=3&mode=trail&max_per_epoch=2")
    assert r.status_code == 200, r.text
    d = r.json()
    assert d["status"] == "ready" and len(d["frames"]) == 3 and d["trail_url"].endswith("/trail.png")
    png = np.asarray(Image.open(cache.CACHE / d["trail_url"].removeprefix("/cache/")).convert("RGB"), float)
    h = png.shape[0]
    # newest epoch (yellow: R and G >> B) is north = upper half; oldest (blue: B >> R) is south = lower half
    top, bottom = png[: h // 2], png[h // 2:]
    yellow = lambda a: (a[..., 0] + a[..., 1] - 2 * a[..., 2]).max()   # noqa: E731
    blue = lambda a: (a[..., 2] - a[..., 0]).max()                       # noqa: E731
    assert yellow(top) > yellow(bottom)
    assert blue(bottom) > blue(top)
    # the oldest epoch is clearly visible, not erased or faded out by the newest (the alpha-blending bug)
    assert blue(bottom) > 60


def test_spectrum_and_lightcurve(fake_irsa):
    r = _wait(f"/api/spectrum?ra={RA}&dec={DEC}&per_epoch=1")
    assert r.status_code == 200, r.text
    s = r.json()
    pts = s["points"]
    assert len(pts) == 6 * 3                                     # 6 detectors x 3 epochs x 1 exposure
    assert [p["wave_um"] for p in pts] == sorted(p["wave_um"] for p in pts)
    assert {p["det"] for p in pts} == set(range(1, 7))
    # the aperture is fixed at RA/Dec (no proper motion for raw coordinates): full flux only at the epoch where
    # the star sits there (epoch index 1), much less at the others
    by_epoch = {}
    for p in pts:
        by_epoch.setdefault(round(p["mjd"]), []).append(p["flux_ujy"])
    assert np.median(by_epoch[61080]) > 5 * np.median(by_epoch[60900])
    n_before = fake_irsa["cutout"]

    lc = c.get(f"/api/lightcurve?ra={RA}&dec={DEC}&per_epoch=1&det=2").json()
    lo, hi = lc["wave_window_um"]
    assert lc["status"] == "ready" and lo < 1.36 < hi
    assert all(lo <= p["wave_um"] <= hi for p in lc["points"])
    assert [p["mjd"] for p in lc["points"]] == sorted(p["mjd"] for p in lc["points"])
    assert fake_irsa["cutout"] == n_before                        # light curve reuses the spectrum, no downloads


def test_lightcurve_bad_window():
    r = c.get(f"/api/lightcurve?ra={RA}&dec={DEC}&wave_min=2.0")
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"
