"""Тесты Exposure.

Офлайн:  uv run pytest -q
С сетью: SPHEREX_LIVE=1 uv run pytest -q   (качает 1 вырезку Барнарда из IRSA)
         PowerShell: $env:SPHEREX_LIVE=1; uv run pytest -q
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from Exposure import Exposure, exposure_from_hdul, load_exposure, query_spherex, cutout_url  # noqa: E402

BARNARD_2025 = SkyCoord(269.4463, 4.768, unit="deg")  # положение Барнарда около 2025.9


def _fake_hdul(n=21):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [BARNARD_2025.ra.deg, BARNARD_2025.dec.deg]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-6.15 / 3600, 6.15 / 3600]
    hdr = w.to_header()
    hdr["DETECTOR"] = 3
    hdr["OBSID"] = "TEST_0001_1"
    hdr["MJD-AVG"] = 60902.5

    img = np.ones((n, n), np.float32)
    img[n // 2, n // 2] = 100.0
    flags = np.zeros((n, n), np.int32)
    flags[0, 0] = (1 << 1) | (1 << 21)   # OVERFLOW + SOURCE
    flags[1, 1] = 1 << 21                # только SOURCE — НЕ плохой
    fh = fits.Header()
    fh["MSKN1"], fh["MSKN21"] = "OVERFLOW", "SOURCE"

    return fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(img, hdr, name="IMAGE"),
        fits.ImageHDU(flags, fh, name="FLAGS"),
        fits.ImageHDU(np.full((n, n), 0.01, np.float32), name="VARIANCE"),
        fits.ImageHDU(np.full((n, n), 0.5, np.float32), name="ZODI"),
    ])


def test_offline_fields_and_masks():
    with pytest.warns(UserWarning, match="WCS-WAVE"):
        e = exposure_from_hdul(_fake_hdul(), source="fake")
    assert isinstance(e, Exposure)
    assert e.detector == 3 and e.mjd == 60902.5 and e.obs_id == "TEST_0001_1"
    x, y = e.sky_to_pixel(BARNARD_2025)
    assert abs(x - 10) < 1e-6 and abs(y - 10) < 1e-6
    assert e.contains(BARNARD_2025, margin_px=3)
    bad = e.bad_mask()
    assert bad[0, 0] and not bad[1, 1]               # SOURCE сам по себе не маскируется
    assert e.flag_names_at(0, 0) == ["OVERFLOW", "SOURCE"]
    assert np.allclose(e.background_subtracted()[5, 5], 0.5)
    wl, _ = e.wavelength_at(x, y)                    # без WCS-WAVE — середина детектора 3
    assert 1.63 < wl < 2.41
    assert 800 < Exposure.mjysr_to_ujy_per_pixel() < 1000


@pytest.mark.skipif(not os.environ.get("SPHEREX_LIVE"), reason="нужна сеть: SPHEREX_LIVE=1")
def test_live_barnard_cutout():
    tbl = query_spherex(BARNARD_2025, detector=6)
    assert len(tbl) > 0
    row = tbl[0]
    e = load_exposure(cutout_url(str(row["access_url"]), BARNARD_2025, 0.05),
                      mjd_hint=float(row["mjd"]))
    assert e.detector == 6 and e.wave_wcs is not None
    assert abs(e.mjd - float(row["mjd"])) < 0.01     # MJD-AVG из заголовка ≈ время из архива
    wl, _ = e.wavelength_at_sky(BARNARD_2025)
    assert 4.42 <= wl <= 5.0
    assert e.contains(BARNARD_2025)
