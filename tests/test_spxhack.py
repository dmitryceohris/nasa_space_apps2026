"""spxhack on synthetic exposures: arrays out, tracking, proper motion with parallax, files. No network."""
import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.time import Time  # noqa: E402

import spxhack as spx  # noqa: E402
from Exposure.exposure import DEFAULT_FLAG_BITS  # noqa: E402
from test_processing import CENTER, fake_exposure, offset  # noqa: E402

RADEC = (CENTER.ra.deg, CENTER.dec.deg)


def mover(pmra=-4000.0, pmdec=15000.0, plx=500.0):
    """A fast nearby star at the test field (Barnard-like, but faster so a few epochs show it)."""
    return spx.Target("fast", "Fast star", CENTER.ra.deg, CENTER.dec.deg, pmra=pmra, pmdec=pmdec,
                      parallax=plx, rv=0.0, ref_epoch=Time(61000, format="mjd"))


def exposures_of(t, epoch_mjds, per=3, flux=50.0, extra=None):
    """per exposures per epoch, the star exactly where spx.predict says (space motion + parallax)."""
    out = []
    for k, m0 in enumerate(epoch_mjds):
        for i in range(per):
            e = fake_exposure(spx.predict(t, m0 + i), m0 + i, flux_mjysr=flux, seed=10 * k + i)
            if extra is not None:
                extra(e)
            out.append(e)
    return out


# ------------------------------------------------------------------ data
def test_target_forms():
    assert spx.target("Barnard's Star").id == "barnard"
    t = spx.target(RADEC)
    assert t.pmra == 0 and t.at(61000).separation(CENTER).arcsec < 1e-6
    assert spx.target(CENTER).ra == pytest.approx(CENTER.ra.deg)
    m = mover()
    assert spx.target(m) is m
    with pytest.raises(TypeError):
        spx.target(42)


def test_exposure_helpers():
    e = fake_exposure(CENTER, 61000)
    e.flags[3, 4] = 1 << DEFAULT_FLAG_BITS["HOT"]
    a = e.to_numpy()
    assert isinstance(a, np.ndarray) and np.isnan(a[3, 4]) and np.isfinite(a[0, 0])
    assert abs(np.nanmedian(a)) < 0.05  # zodi removed
    assert e.det == 3 and 1.63 < e.wave(10, 10) < 2.41


def test_epochs_from_exposures():
    exps = [fake_exposure(CENTER, m) for m in (100, 101, 300)]
    assert spx.epochs(exps) == [[0, 1], [2]]
    assert spx.epochs(np.array([100.0, 101.0, 300.0])) == [[0, 1], [2]]


def test_search_and_fetch_offline(monkeypatch):
    tbl = Table({"obs_id": ["a", "b"], "detector": [3, 3], "mjd": [100.0, 300.0], "em_min": [1.9e-6, 2.0e-6],
                 "em_max": [1.95e-6, 2.05e-6], "access_url": ["u1", "u2"]})
    monkeypatch.setattr("spxhack.data.query_spherex", lambda c, detector=None: tbl)
    rows = spx.search("barnard", det=3)
    assert [r["obs_id"] for r in rows] == ["a", "b"] and rows[0]["wave_min_um"] == pytest.approx(1.9)
    got = {}
    monkeypatch.setattr("spxhack.data.download_rows",
                        lambda rows, c, size, progress=None: got.setdefault("n", len(rows)) and [])
    spx.fetch("barnard", det=3, fov=3)
    assert got["n"] == 2


# ------------------------------------------------------------------ stack and arrays
def test_stack_arrays_and_coordinates():
    star = offset(CENTER, 5, -4)
    exps = [fake_exposure(star, 60900 + 200 * k + i, seed=10 * k + i) for k in range(3) for i in range(3)]
    st = spx.stack(exps, RADEC, fov=3)
    assert len(st) == 3 and st.cube.shape == (3, 29, 29) == np.asarray(st).shape
    assert st.to_numpy() is st.cube and st.variance.shape == st.flags.shape == st.cube.shape
    assert st.mjd.shape == (3,) and len(st.dates) == 3 and (st.n_exp == 3).all()
    x, y = st.pixel(star)
    ra, dec = st.radec(x, y)
    assert ra == pytest.approx(star.ra.deg) and dec == pytest.approx(star.dec.deg)
    assert st.pixel_scale == pytest.approx(6.15)
    d, info = spx.diff(st, info=True)
    assert d.shape == st.cube.shape and len(info["noise_scale"]) == 3 and info["method"]["name"]
    assert spx.trail(st).shape == (29, 29, 3)


def test_single_epoch_has_no_diff():
    st = spx.stack([fake_exposure(CENTER, 61000)], RADEC)
    with pytest.raises(ValueError):
        spx.diff(st)
    with pytest.raises(ValueError):
        spx.trail(st)
    with pytest.raises(ValueError):
        spx.stack([], RADEC)


def test_centroid_on_plain_array():
    yy, xx = np.indices((30, 30))
    img = np.random.default_rng(0).normal(0, 1, (30, 30)) + 40 * np.exp(-((xx - 14.3) ** 2 + (yy - 15.6) ** 2) / 4)
    m = spx.centroid(img.tolist())  # lists work too
    assert m["status"] == "ok" and abs(m["xy"][0] - 14.3) < 0.3 and abs(m["xy"][1] - 15.6) < 0.3


# ------------------------------------------------------------------ track and proper motion
EPOCHS5 = [60900, 61000, 61090, 61180, 61270]  # ~1 year, sampling the parallax ellipse


def test_track_follows_prediction():
    t = mover()
    st = spx.stack(exposures_of(t, EPOCHS5), t, fov=4)
    tr = spx.track(st)
    assert len(tr) == 5 and all(r["association"] == "unique" for r in tr)
    assert all(r["offset_px"] < 0.2 for r in tr)
    a = tr.to_numpy()
    assert a.shape == (5, len(spx.Track.columns)) and np.isfinite(a).all()
    assert np.asarray(tr).shape == a.shape


def test_fit_motion_recovers_catalogue_pm_with_parallax():
    t = mover()
    st = spx.stack(exposures_of(t, EPOCHS5), t, fov=4)
    m = spx.fit_motion(spx.track(st))
    assert m.parallax_mode == "catalog" and m.n_epochs == 5
    assert m.to_numpy() == pytest.approx([-4.0, 15.0], abs=0.3)
    assert np.all(np.abs(m.difference_sigma) < 5)
    assert m.to_dict()["pm"] == pytest.approx(list(m.pm))
    # ignoring a 0.5" parallax over 5 epochs biases the fit more than fixing it
    m0 = spx.fit_motion(spx.track(st), parallax="none")
    assert np.hypot(*(m0.pm - [-4, 15])) > np.hypot(*(m.pm - [-4, 15]))
    mf = spx.fit_motion(spx.track(st), parallax="fit")
    assert mf.parallax_mas == pytest.approx(500, abs=200)


def test_fit_motion_needs_two_epochs_and_valid_mode():
    t = mover()
    st = spx.stack(exposures_of(t, [60900]), t)
    with pytest.raises(ValueError):
        spx.fit_motion(spx.track(st))
    st = spx.stack(exposures_of(t, [60900, 61100]), t)
    with pytest.raises(ValueError):
        spx.fit_motion(spx.track(st), parallax="fit")  # needs 4 epochs
    with pytest.raises(ValueError):
        spx.fit_motion(spx.track(st), parallax="maybe")
    m = spx.fit_motion(spx.track(st))
    assert m.chi2_red is None and m.warnings  # 2 epochs: no degrees of freedom, said so


def _add_star(where, amp, sigma_px):
    def add(e):
        x0, y0 = e.wcs.world_to_pixel(where)
        yy, xx = np.indices(e.shape)
        e.image += (amp * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma_px**2))).astype(np.float32)
    return add


def test_ambiguous_frame_resolved_by_prediction():
    """A bigger, separate star off-centre makes every frame ambiguous; the target is picked by its prediction."""
    t = mover(pmra=0, pmdec=8000, plx=0)
    exps = exposures_of(t, [60900, 61100, 61300], extra=_add_star(offset(CENTER, 75, 60), 30, 1.8))
    tr = spx.track(spx.stack(exps, t, fov=5))
    assert all(r["status"] == "ambiguous" for r in tr)
    assert all(r["association"] == "nearest_to_prediction" and r["offset_px"] < 0.3 for r in tr)
    assert spx.fit_motion(tr).pm[1] == pytest.approx(8.0, abs=0.5)


def test_blend_with_bright_neighbour_measured_at_prediction():
    """The neighbour's wings touch the target: one merged blob centred between them. Windowed centroid at the
    prediction recovers the target and says it was blended."""
    t = mover(pmra=0, pmdec=8000, plx=0)
    exps = exposures_of(t, [60900, 61100, 61300], extra=_add_star(offset(CENTER, 40, 40), 60, 3.0))
    st = spx.stack(exps, t, fov=4)
    assert spx.centroid(st.cube[0])["candidates"][0]["area"] > 100  # really merged
    tr = spx.track(st)
    for r in tr:
        assert r["association"] == "windowed_at_prediction" and "blended" in r["flags"]
        assert r["offset_px"] < 0.3 and r["quality"] == "fair"
    assert spx.fit_motion(tr).pm[1] == pytest.approx(8.0, abs=0.3)  # 8.9 without the local plane background


def test_undetected_target_forced_only_when_significant():
    t = mover(pmra=0, pmdec=8000, plx=0)
    faint = spx.stack(exposures_of(t, [60900, 61100], flux=0.0), t, fov=3)  # nothing there
    assert all(r["xy"] is None for r in spx.track(faint))


def test_far_detection_is_flagged_and_not_fitted():
    """The catalogue says the target is here, but the only object is 6 px away: not the target."""
    t = mover(pmra=0, pmdec=0, plx=100)
    other = offset(CENTER, -37, 0)
    exps = [fake_exposure(other, 60900 + 200 * k + i, seed=10 * k + i) for k in range(3) for i in range(3)]
    tr = spx.track(spx.stack(exps, t, fov=4))
    assert all("far_from_prediction" in r["flags"] for r in tr)
    assert tr.good() == []
    with pytest.raises(ValueError):
        spx.fit_motion(tr)


def test_parallax_offset_amplitude():
    mjd = np.linspace(61000, 61365, 50)
    d = spx.parallax_offset(270.0, 66.56, mjd, 1000.0)  # ecliptic pole: a circle of radius = parallax
    r = np.hypot(d[:, 0], d[:, 1])
    assert r.min() > 0.95 and r.max() < 1.05
    assert spx.parallax_offset(10.0, 0.0, 61000.0, 0.0) == pytest.approx([0, 0])


# ------------------------------------------------------------------ photometry
def test_spectrum_and_lightcurve_to_numpy():
    exps = [fake_exposure(CENTER, 60900 + 100 * i, flux_mjysr=40, sigma_px=1.2, noise=0.01, seed=i)
            for i in range(4)]
    sp = spx.spectrum(exps, RADEC)
    a = sp.to_numpy()
    assert a.shape == (4, 3) and np.all(a[:, 1] > 0) and sp.ok.all()
    assert np.asarray(sp).shape == (4, 3) and sp.method["name"] == "forced_aperture_photometry"
    lc = spx.lightcurve(sp, 1.9, 2.2)
    b = lc.to_numpy()
    assert b.shape == (4, 3) and np.all(np.diff(b[:, 0]) > 0)


# ------------------------------------------------------------------ output
def test_to_fits_round_trip(tmp_path):
    t = mover()
    st = spx.stack(exposures_of(t, [60900, 61100, 61300]), t, fov=3)
    p = spx.to_fits(st, tmp_path / "c.fits")
    with fits.open(p) as h:
        assert h[0].data.shape == st.cube.shape and h[0].header["BUNIT"] == "MJy/sr"
        assert np.allclose(h[0].data, st.cube.astype(np.float32), equal_nan=True)
        assert [x.name for x in h[1:]] == ["VARIANCE", "FLAGS", "EPOCHS"] and len(h["EPOCHS"].data) == 3
        assert h[0].header["CTYPE1"] == "RA---TAN"
    p2 = spx.to_fits(st, tmp_path / "d.fits", what="diff")
    assert fits.getheader(p2)["BUNIT"] == "sigma"
    with pytest.raises(OSError):
        spx.to_fits(st, p)  # no silent overwrite


def test_save_and_show(tmp_path):
    t = mover()
    st = spx.stack(exposures_of(t, [60900, 61100, 61300]), t, fov=3)
    res = spx.save(st, tmp_path / "out")
    assert len(res["frames"]) == 3 and (tmp_path / "out" / "result.json").exists()
    from matplotlib.animation import FuncAnimation

    for mode in ("coadd", "diff"):
        anim = spx.show(st, mode=mode, track=spx.track(st))
        assert isinstance(anim, FuncAnimation)
        assert "<img" in anim.to_jshtml() or "frames" in anim.to_jshtml()
    with pytest.raises(ValueError):
        spx.show(st, mode="rgb")
