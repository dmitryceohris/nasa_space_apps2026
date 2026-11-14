"""Тесты обработки на синтетических снимках (без сети)."""
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
from PIL import Image

from Exposure import Exposure
from processing import (coadd, difference, find_preset, forced_photometry, group_epochs,
                        load_targets, make_grid, render_view, select_per_epoch, trail)
from processing.frames import grid_header_json
from processing.render import fits_to_png_xy

CENTER = SkyCoord(150.0, 20.0, unit="deg")
PIX = 6.15


def fake_exposure(star: SkyCoord | None, mjd: float, flux_mjysr: float = 50.0, n: int = 41,
                  sigma_px: float = 1.0, noise: float = 0.05, seed: int = 0, shift_px=(0.3, -0.2)) -> Exposure:
    """Снимок с одной гауссовой звездой; WCS слегка сдвинут относительно общей сетки."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    w.wcs.crpix = [(n + 1) / 2 + shift_px[0], (n + 1) / 2 + shift_px[1]]
    w.wcs.cdelt = [-PIX / 3600, PIX / 3600]
    rng = np.random.default_rng(seed)
    img = rng.normal(0, noise, (n, n)).astype(np.float32) + 1.0  # +1 = «зодиакальный фон»
    if star is not None:
        x0, y0 = w.world_to_pixel(star)
        yy, xx = np.indices((n, n))
        img += (flux_mjysr * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma_px ** 2))).astype(np.float32)
    return Exposure(image=img, variance=np.full((n, n), noise ** 2, np.float32),
                    flags=np.zeros((n, n), np.int32), zodi=np.ones((n, n), np.float32),
                    wcs=w, wave_wcs=None, detector=3, mjd=mjd, obs_id=f"T{mjd}")


def offset(c: SkyCoord, east_arcsec: float, north_arcsec: float) -> SkyCoord:
    return c.spherical_offsets_by(east_arcsec * u.arcsec, north_arcsec * u.arcsec)


# ------------------------------------------------------------------ эпохи
def test_group_epochs_by_gap_and_span():
    mjds = [100, 101, 103, 300, 302, 500]
    assert group_epochs(mjds) == [[0, 1, 2], [3, 4], [5]]
    long = list(range(0, 60, 2))  # 60 дней непрерывно (полюс эклиптики) -> несколько кадров
    assert len(group_epochs(long, max_span_days=20)) == 3


def test_select_per_epoch_even():
    mjds = list(range(10)) + [200]
    sel = select_per_epoch(mjds, 3)
    assert sel[0] == [0, 4, 9] and sel[1] == [10]


# ------------------------------------------------------------------ сетка и кадры
def test_grid_center_and_orientation():
    wcs, shape = make_grid(CENTER, 3.0)
    assert shape[0] % 2 == 1 and shape == (shape[0], shape[0])
    x, y = wcs.world_to_pixel(CENTER)
    assert abs(x - (shape[1] - 1) / 2) < 1e-6 and abs(y - (shape[0] - 1) / 2) < 1e-6
    _, yn = wcs.world_to_pixel(offset(CENTER, 0, 30))
    xe, _ = wcs.world_to_pixel(offset(CENTER, 30, 0))
    assert yn > y and xe < x  # север — вверх по оси y FITS, восток — влево
    assert grid_header_json(wcs, shape)["pixel_scale_arcsec"] == pytest.approx(PIX)


def test_moving_star_gives_plus_minus_pair():
    old, new = offset(CENTER, 0, -15), offset(CENTER, 0, 15)  # сдвиг на ~5 пикселей к северу
    e1 = [fake_exposure(old, 60900 + i, seed=i) for i in range(3)]
    e2 = [fake_exposure(new, 61100 + i, seed=10 + i) for i in range(3)]
    e3 = [fake_exposure(new, 61300 + i, seed=20 + i) for i in range(3)]
    wcs, shape = make_grid(CENTER, 3.0)
    frames = [coadd(e, wcs, shape, CENTER) for e in (e1, e2, e3)]
    assert all(f.n_exp == 3 for f in frames)
    diffs, _ = difference(frames)
    xo, yo = (int(round(float(v))) for v in wcs.world_to_pixel(old))
    xn, yn = (int(round(float(v))) for v in wcs.world_to_pixel(new))
    assert diffs[0][yo, xo] > 10       # в первой эпохе звезда там, где шаблона нет: «+»
    assert diffs[0][yn, xn] < -10      # и нет там, где она в шаблоне: «−»
    rgb = trail(diffs)
    assert rgb.shape == (*shape, 3) and 0 <= rgb.min() and rgb.max() <= 1


def test_coadd_reduces_noise_and_removes_background():
    wcs, shape = make_grid(CENTER, 3.0)
    one = coadd([fake_exposure(None, 1, seed=1)], wcs, shape, CENTER)
    four = coadd([fake_exposure(None, i, seed=i) for i in range(4)], wcs, shape, CENTER)
    assert abs(np.nanmedian(four.image)) < 0.01                     # фон вычтен
    assert np.nanstd(four.image) < 0.7 * np.nanstd(one.image)       # шум ~ 1/√N


# ------------------------------------------------------------------ рендер
def test_png_has_north_up(tmp_path):
    star = offset(CENTER, 0, 60)  # звезда севернее центра -> в верхней половине PNG
    exps = [fake_exposure(star, 60900, seed=1), fake_exposure(star, 61100, seed=2)]
    wcs, shape = make_grid(CENTER, 3.0)
    frames = [coadd([e], wcs, shape, CENTER) for e in exps]
    res = render_view(frames, wcs, shape, tmp_path, "/c")
    png = np.asarray(Image.open(tmp_path / "f0.png"))
    row = np.unravel_index(np.argmax(png), png.shape)[0]
    assert row < shape[0] // 2
    x, y = wcs.world_to_pixel(star)
    assert fits_to_png_xy(x, y, shape[0])[1] == pytest.approx(row, abs=1)
    assert len(res["frames"]) == 2 and res["trail_url"] == "/c/trail.png"
    for name in ("f0.png", "f0_diff.png", "f0_flags.png", "f1.png", "trail.png", "movie.gif"):
        assert (tmp_path / name).exists()


def test_render_view_empty(tmp_path):
    wcs, shape = make_grid(CENTER, 2.0)
    res = render_view([], wcs, shape, tmp_path, "/c")
    assert res["frames"] == [] and res["warnings"]


# ------------------------------------------------------------------ фотометрия
def test_forced_photometry_recovers_flux():
    sigma, amp = 1.2, 40.0
    e = fake_exposure(CENTER, 61000, flux_mjysr=amp, sigma_px=sigma, noise=0.01)
    p = forced_photometry(e, CENTER, r_ap=4.0, r_in=7, r_out=12)
    true_ujy = amp * 2 * np.pi * sigma ** 2 * Exposure.mjysr_to_ujy_per_pixel()
    assert p["ok"] and p["flux_ujy"] == pytest.approx(true_ujy, rel=0.03)
    assert 1.63 < p["wave_um"] < 2.41 and p["det"] == 3


# ------------------------------------------------------------------ цели
def test_targets_and_motion():
    t = load_targets()
    assert set(t) == {"barnard", "epsindib", "wise0855", "rwaur", "rdra"}
    b = t["barnard"]
    c1, c2 = b.coord_at(Time(2025.5, format="jyear")), b.coord_at(Time(2026.5, format="jyear"))
    assert c1.separation(c2).arcsec == pytest.approx(10.39, abs=0.05)   # ~10.4″ за год
    w = t["wise0855"]
    assert w.ref_epoch.mjd == pytest.approx(57633.39)
    assert find_preset("Barnard's Star", t) is b and find_preset("WISE 0855-0714", t) is w
