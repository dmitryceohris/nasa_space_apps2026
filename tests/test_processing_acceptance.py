"""B1 acceptance through the user path GET /api/view (synthetic exposures replace only the IRSA download).

Each test is one criterion of DILAYLA_PLAN.md B1 / docs/reviews/processing-acceptance.md."""
import os

os.environ["PREWARM"] = "0"

import time  # noqa: E402

import numpy as np  # noqa: E402
import pytest  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402
from PIL import Image  # noqa: E402

from api import cache, irsa, jobs, resolve  # noqa: E402
from api.main import app  # noqa: E402
from Exposure.exposure import DEFAULT_FLAG_BITS  # noqa: E402
from processing.render import NODATA_RGB  # noqa: E402
from test_processing import CENTER, PIX, fake_exposure, offset  # noqa: E402

SOURCE, HOT = 1 << DEFAULT_FLAG_BITS["SOURCE"], 1 << DEFAULT_FLAG_BITS["HOT"]


@pytest.fixture
def serve(monkeypatch, tmp_path):
    """serve(exposures, mode, fov=3, motion=None) -> /api/view JSON for exactly these exposures."""
    monkeypatch.setattr(cache, "CACHE", tmp_path)
    client = TestClient(app)

    def run(exps, mode="single", fov=3.0, max_per_epoch=6, sky=None, failing=()):
        meta = [{"uri": f"u{i}", "mjd": e.mjd, "lam": [2.0, 2.1]} for i, e in enumerate(exps)]
        by_uri = {f"u{i}": e for i, e in enumerate(exps)}
        monkeypatch.setattr(irsa, "exposures", lambda ra, dec, det: meta)

        def cutout(uri, *a):
            if uri in failing:
                raise OSError("controlled download failure")
            return uri
        monkeypatch.setattr(irsa, "cutout", cutout)
        monkeypatch.setattr(jobs, "load_exposure", lambda path, mjd_hint=None, obs_id_hint="": by_uri[path])
        if sky is not None:  # a moving target: prediction per date
            monkeypatch.setattr(resolve, "for_params", lambda p: {
                "center": (CENTER.ra.deg, CENTER.dec.deg), "label": "controlled target",
                "motion": {"name": "controlled target", "source": "fixture", "pmra_masyr": 0, "pmdec_masyr": 1}})
            monkeypatch.setattr(jobs, "_sky", lambda r, mjd: sky(mjd))
        q = {"ra": CENTER.ra.deg, "dec": CENTER.dec.deg, "fov": fov, "det": 3, "mode": mode,
             "max_per_epoch": max_per_epoch}
        for _ in range(400):
            r = client.get("/api/view", params=q)
            if r.status_code != 202:
                break
            time.sleep(0.01)
        assert r.status_code == 200, r.text
        res = r.json()
        res["_dir"] = tmp_path
        return res
    return run


def png(res, i, tmp):
    return np.asarray(Image.open(tmp / res["frames"][i]["image_url"].removeprefix("/cache/")).convert("RGB"))


def png_xy(res, coord):
    """Expected PNG pixel of a sky position on the response grid."""
    from astropy.wcs import WCS
    g = res["grid"]
    w = WCS(naxis=2)
    w.wcs.ctype = [g["wcs"]["CTYPE1"], g["wcs"]["CTYPE2"]]
    w.wcs.crval = [g["wcs"]["CRVAL1"], g["wcs"]["CRVAL2"]]
    w.wcs.crpix = [g["wcs"]["CRPIX1"], g["wcs"]["CRPIX2"]]
    w.wcs.cdelt = [g["wcs"]["CDELT1"], g["wcs"]["CDELT2"]]
    x, y = w.world_to_pixel(coord)
    return float(x), float(g["height_px"] - 1 - y)


def dist(a, b):
    return float(np.hypot(a[0] - b[0], a[1] - b[1]))


# ------------------------------------------------------------------ 1. masks before interpolation, SOURCE, coverage
@pytest.mark.parametrize("flagged", [True, False])
def test_single_bad_native_pixel_never_becomes_the_target(serve, flagged):
    """One bad native pixel in an empty field. Flagged (HOT) or not (missed cosmic ray), it must not become a
    source: bilinear reprojection used to spread an unflagged one into a S/N ~280 compact 'star'."""
    exps = []
    for k, mjd in enumerate([61000.0, 61200.0]):
        e = fake_exposure(None, mjd, seed=k)
        e.image[24, 27] += 50.0  # 1000 sigma, one native pixel
        if flagged:
            e.flags[24, 27] |= HOT
        exps.append(e)
    res = serve(exps, "single")
    assert len(res["frames"]) == 2
    for f in res["frames"]:
        assert f["measured_px"] is None and f["measured_candidates"] == []


def test_source_flag_survives_the_native_spike_heuristic_through_http(serve, monkeypatch):
    import reproject
    original = reproject.reproject_interp
    seen = []

    def project(input_data, *args, **kwargs):
        arr, _ = input_data
        if kwargs.get("order") == "bilinear":
            assert np.isfinite(arr[24, 27]), "SOURCE alone must not become a missing pixel"
            seen.append(True)
        return original(input_data, *args, **kwargs)

    monkeypatch.setattr(reproject, "reproject_interp", project)
    exps = []
    for seed, mjd in enumerate([61000., 61200.]):
        exp = fake_exposure(None, mjd, seed=seed)
        exp.image[24, 27] += 50
        exp.flags[24, 27] = SOURCE
        exps.append(exp)
    result = serve(exps)
    assert len(result["frames"]) == 2 and len(seen) == 4


def test_source_is_kept_and_no_data_has_its_own_colour(serve):
    """SOURCE marks the object, it is never masked; pixels without coverage are not painted as dark sky."""
    star = offset(CENTER, 0, 0)
    exps = []
    for k, mjd in enumerate([61000.0, 61200.0]):
        e = fake_exposure(star, mjd, seed=k, shift_px=(9.0, 0.3))  # the exposure covers only part of the grid
        e.flags[18:23, 18:23] |= SOURCE
        x0, y0 = (int(round(v)) for v in e.sky_to_pixel(star))
        e.flags[y0 - 2:y0 + 3, x0 - 2:x0 + 3] |= SOURCE
        exps.append(e)
    res = serve(exps, "single")
    f = res["frames"][0]
    assert f["quality"]["coverage_fraction"] < 1
    assert f["measured_px"] is not None and dist(f["measured_px"], png_xy(res, star)) < 0.3
    rgb = png(res, 0, res["_dir"])
    nodata = np.all(rgb == NODATA_RGB, axis=2)
    grey = (rgb[..., 0] == rgb[..., 1]) & (rgb[..., 1] == rgb[..., 2])
    assert nodata.any() and grey.any()
    assert np.all(nodata | grey)  # every pixel is either sky/source (grey) or explicitly "no data"


@pytest.mark.parametrize("mode", ["coadd", "diff", "trail"])
def test_missing_pixels_keep_their_colour_in_each_view_mode(serve, mode):
    exps = [fake_exposure(CENTER, mjd, seed=i, shift_px=(9.0, .3))
            for i, mjd in enumerate([61000., 61200.])]
    result = serve(exps, mode)
    rgb = png(result, 0, result["_dir"])
    assert np.all(rgb == NODATA_RGB, axis=2).any()
    if mode in ("diff", "trail"):
        path = result["_dir"] / result["trail_url"].removeprefix("/cache/")
        rgb = np.asarray(Image.open(path).convert("RGB"))
        assert np.all(rgb == NODATA_RGB, axis=2).any()


# ------------------------------------------------------------------ 2. background and NaN edge
def test_stationary_star_with_wrong_zodi_models_shows_no_change(serve):
    """Zodi model off by a different constant and gradient in each epoch: the star must not move or
    leave a residual pair; the coadd centre stays put."""
    star = offset(CENTER, 3, -2)
    errs = [(0.3, 0.0), (-0.25, 0.004), (0.1, -0.006)]
    exps = []
    for k, (c, g) in enumerate(errs):
        for i in range(2):
            e = fake_exposure(star, 61000.0 + 200 * k + i, flux_mjysr=5.0, seed=10 * k + i)
            yy, xx = np.indices(e.shape)
            e.zodi = e.zodi - c - g * (xx - 20)  # the model misses c + gradient
            exps.append(e)
    co = serve(exps, "coadd")
    true = png_xy(co, star)
    xs = [f["measured_px"] for f in co["frames"]]
    assert all(x is not None and dist(x, true) < 0.15 for x in xs)
    assert all(abs(o - c) < 0.05 for f, (c, g) in zip(co["frames"], errs) for o in f["bg_offsets"])
    diff = serve(exps, "diff")
    for f in diff["frames"]:
        plus = [c for c in f["measured_candidates"] if c["sign"] == 1 and dist(c["xy"], true) < 3]
        assert f["measured_px"] is None and plus == []


def test_nan_edge_next_to_target_does_not_pull_the_centre(serve):
    star = offset(CENTER, 0, 0)
    exps = []
    for k, mjd in enumerate([61000.0, 61200.0]):
        e = fake_exposure(star, mjd, seed=k)
        x0, _ = e.sky_to_pixel(star)
        e.image[:, int(x0) + 4:] = np.nan  # a no-data strip 4 px from the star (e.g. detector edge)
        exps.append(e)
    res = serve(exps, "single")
    true = png_xy(res, star)
    for f in res["frames"]:
        c = f["measured_candidates"]
        assert len(c) == 1 and dist(c[0]["xy"], true) < 0.3
        # the old centroid pulled the centre several px into the strip; now it is at the star and, because the
        # star touches missing data, it is reported as a candidate only
        assert "touches_nodata" in c[0]["flags"] and f["measured_px"] is None


# ------------------------------------------------------------------ 3. epochs: duration and real composition
def test_long_series_is_split_and_times_are_those_of_used_exposures(serve):
    exps = [fake_exposure(CENTER, 61000.0 + 5 * i, seed=i) for i in range(12)]  # 55 days, every 5 days
    res = serve(exps, "coadd", failing={"u1"})  # one download fails
    assert len(res["frames"]) >= 3
    used = set()
    for f in res["frames"]:
        mjds = [float(o[1:]) for o in f["obs_ids"]]
        assert f["mjd_span"][1] - f["mjd_span"][0] <= 20
        assert f["mjd"] == pytest.approx(np.mean(mjds), abs=1e-3) and f["n_exp"] == len(mjds)
        used |= set(f["obs_ids"])
    assert "T61005.0" not in used and len(used) == 11
    assert res["exposure_counts"]["failed"] == 1


# ------------------------------------------------------------------ 4. diff: separate + and -
def test_sub_psf_motion_in_diff_is_not_reported_as_a_position(serve):
    """8" = 1.3 px: '+' and '-' touch. Measured on '+' alone the centre is still pulled ~0.9 px away from the
    '-' (2/3 of the motion), so the frame gives candidates, not a position; the coadd frames give positions."""
    old, new = offset(CENTER, 0, -4), offset(CENTER, 0, 4)
    exps = [fake_exposure(old, 61000.0 + i, seed=i) for i in range(2)] +            [fake_exposure(new, 61200.0 + i, seed=10 + i) for i in range(2)]
    res = serve(exps, "diff")
    for f in res["frames"]:
        assert f["measured_px"] is None
        assert "overlaps_negative_residual" in f["measurement"]["flags"]
        assert {c["sign"] for c in f["measured_candidates"]} == {1, -1}  # both parts listed separately
    co = serve(exps, "coadd")
    assert dist(co["frames"][0]["measured_px"], png_xy(co, old)) < 0.2
    assert dist(co["frames"][1]["measured_px"], png_xy(co, new)) < 0.2


def test_resolved_motion_is_measured_on_the_positive_part(serve):
    old, new = offset(CENTER, 0, -18), offset(CENTER, 0, 18)  # 36" = 5.9 px apart
    exps = [fake_exposure(old, 61000.0 + i, seed=i) for i in range(2)] +            [fake_exposure(new, 61200.0 + i, seed=10 + i) for i in range(2)]
    res = serve(exps, "diff")
    p_old, p_new = png_xy(res, old), png_xy(res, new)
    f0, f1 = res["frames"]
    assert dist(f0["measured_px"], p_old) < 0.4 and dist(f1["measured_px"], p_new) < 0.4
    neg = [c for c in f0["measured_candidates"] if c["sign"] == -1]
    assert neg and dist(neg[0]["xy"], p_new) < 0.4  # not merged into one centre between the two
    mid = ((p_old[0] + p_new[0]) / 2, (p_old[1] + p_new[1]) / 2)
    assert dist(f0["measured_px"], mid) > 2.5
    assert res["processing"]["difference"]["units"] == "robust_residual_scale"


# ------------------------------------------------------------------ 7. spectral channels
@pytest.mark.parametrize("waves,flagged", [((2.00, 2.06), True), ((2.000, 2.010), False)])
def test_different_spectral_channels_are_not_hidden(serve, waves, flagged):
    """3% apart is below the old 5% rule but two channels apart at R~41: the user must be told."""
    exps = []
    for k, w in enumerate(waves):
        e = fake_exposure(CENTER, 61000.0 + 200 * k, seed=k)
        e.wavelength_at = (lambda w: lambda x, y: (w, 0.049))(w)
        exps.append(e)
    res = serve(exps, "single")
    assert res["spectral_sampling"]["different_channels"] is flagged
    assert any("разных спектральных каналах" in w for w in res["warnings"]) is flagged


# ------------------------------------------------------------------ 8. measured_px only for the target
def test_detection_far_from_the_prediction_is_not_the_target(serve):
    """WISE 0855 at 4.42 um (public-final): the target was not detected, the only object 17 px away was
    returned as measured_px. A centre that disagrees with the prediction stays a candidate."""
    other = offset(CENTER, 0, 60)  # 10 px north of the target
    exps = [fake_exposure(other, 61000.0 + 200 * k, seed=k) for k in range(2)]
    res = serve(exps, "single", sky=lambda mjd: CENTER)
    for f in res["frames"]:
        assert f["predicted_px"] is not None and f["measured_px"] is None
        assert f["measurement"]["target_association"] == "far_from_prediction"
        assert f["measurement"]["offset_from_prediction_px"] > 1 and f["measured_candidates"]


def test_saturated_centre_1_5_px_off_is_not_the_target(serve):
    """Replay of Barnard 2026W19_1A_0046_1: a 'fair' centre 1.5 px from the prediction. Registration is good to
    0.4 px, so 1.5 px is not a measurement of the target."""
    star = offset(CENTER, 0, -9)  # 1.46 px south of the prediction
    exps = [fake_exposure(star, 61000.0 + 200 * k, seed=k) for k in range(2)]
    res = serve(exps, "single", sky=lambda mjd: CENTER)
    for f in res["frames"]:
        assert f["measured_px"] is None and f["measurement"]["target_association"] == "far_from_prediction"


def test_target_at_prediction_is_measured(serve):
    exps = [fake_exposure(CENTER, 61000.0 + 200 * k, seed=k) for k in range(2)]
    res = serve(exps, "single", sky=lambda mjd: CENTER)
    for f in res["frames"]:
        assert f["measured_px"] is not None
        assert f["measurement"]["target_association"] == "consistent_with_prediction"


# ------------------------------------------------------------------ 6. photometry errors and coverage
def test_photometry_error_includes_background_and_ring_coverage_counts():
    from processing.photometry import forced_photometry

    e = fake_exposure(CENTER, 61000.0, flux_mjysr=20, sigma_px=1.2, noise=0.05, seed=1)
    p = forced_photometry(e, CENTER)
    k = e.mjysr_to_ujy_per_pixel()
    x0, y0 = e.sky_to_pixel(CENTER)
    yy, xx = np.indices(e.shape)
    ap = np.hypot(xx - x0, yy - y0) <= 2.0
    pixel_only = np.sqrt(e.variance[ap].sum()) * k
    assert p["ok"] and p["err_ujy"] > 1.05 * pixel_only  # + error of the ring median

    ring = (np.hypot(xx - x0, yy - y0) >= 4) & (np.hypot(xx - x0, yy - y0) <= 7)
    e.flags[ring & (xx < x0 + 3)] |= HOT  # most of the ring masked
    p = forced_photometry(e, CENTER)
    assert not p["ok"] and "background_annulus_mostly_masked" in p["quality_reasons"]


def test_photometry_masks_an_unflagged_spike_in_the_aperture():
    from processing.photometry import forced_photometry

    e = fake_exposure(CENTER, 61000.0, flux_mjysr=2, sigma_px=1.2, noise=0.05, seed=2)
    x0, y0 = (int(round(v)) for v in e.sky_to_pixel(CENTER))
    clean = forced_photometry(e, CENTER)["flux_ujy"]
    e.image[y0 + 1, x0 + 1] += 40.0  # cosmic ray next to the centre, not in FLAGS... but the star is there
    # next to a real star the neighbours are bright, so the spike test keeps the pixel: the result is honest
    # only when the spike is isolated. Isolated spike in the aperture of an empty position:
    e2 = fake_exposure(None, 61000.0, noise=0.05, seed=3)
    e2.image[y0 + 1, x0 + 1] += 40.0
    p2 = forced_photometry(e2, CENTER)
    assert not p2["ok"] and "masked_aperture_pixels" in p2["quality_reasons"]
    assert abs(p2["flux_ujy"]) < 0.1 * clean


def test_lightcurve_mixing_channels_is_flagged():
    from processing.pipeline import lightcurve_from_spectrum

    pts = [{"wave_um": w, "bw_um": 0.05, "flux_ujy": 10.0, "err_ujy": 0.1, "ok": True, "mjd": 61000 + i}
           for i, w in enumerate([2.00, 2.02, 2.09])]
    lc = lightcurve_from_spectrum({"points": pts}, 1.95, 2.10)
    assert lc["spectral_sampling"]["different_channels"] and lc["lightcurve_warnings"]
    same = lightcurve_from_spectrum({"points": pts[:2]}, 1.95, 2.10)
    assert not same["spectral_sampling"]["different_channels"]


def test_pixel_scale_constant():
    assert PIX == pytest.approx(6.15)


# ------------------------------------------------------------------ code review follow-ups
def test_prediction_outside_the_field_gives_no_measured_position(serve):
    """The model puts the target outside the frame: the star at the centre is not the target."""
    exps = [fake_exposure(CENTER, 61000.0 + 200 * k, seed=k) for k in range(2)]
    res = serve(exps, "single", sky=lambda mjd: SkyCoord(151, 20, unit="deg"))
    for f in res["frames"]:
        assert f["prediction_status"] == "outside_field" and f["measured_px"] is None
        assert f["measurement"]["target_association"] == "prediction_outside_field"
        assert f["measured_candidates"]  # still listed as a candidate


def test_trail_keeps_the_no_data_colour(serve):
    exps = []
    for k, mjd in enumerate([61000.0, 61200.0]):
        exps.append(fake_exposure(CENTER, mjd, seed=k, shift_px=(9.0, 0.3)))  # partial coverage
    res = serve(exps, "trail")
    rgb = np.asarray(Image.open(res["_dir"] / res["trail_url"].removeprefix("/cache/")).convert("RGB"))
    nodata = np.all(rgb == NODATA_RGB, axis=2)
    assert nodata.any() and not nodata.all()


def test_dipole_overlap_found_even_when_the_minus_part_is_not_listed():
    """A big '-' residual elsewhere is 'largest', another '-' under the centre is 'central': the target's own
    '-' part next to its '+' is neither, and must still mark the '+' centre as biased."""
    from processing.centroid import measure

    yy, xx = np.indices((61, 61))
    g = lambda x, y, s, a: a * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s ** 2))  # noqa: E731
    d = np.random.default_rng(0).normal(0, 1, (61, 61))
    d += g(45, 30, 1.2, 40) - g(45, 32, 1.2, 40)  # the moving target, off-centre: + and - touch
    d -= g(10, 10, 4.0, 40)  # a big negative residual far away -> 'largest'
    d -= g(30, 30, 1.0, 20)  # a small negative residual at the centre -> 'central'
    m = measure(d, signed=True)
    neg = [c for c in m["candidates"] if c["sign"] == -1]
    assert all(abs(c["xy"][0] - 45) > 3 for c in neg)  # the target's '-' is not among the listed ones
    assert "overlaps_negative_residual" in m["flags"] and m["quality"] == "poor"


@pytest.mark.parametrize("status", [400, 403, 404])
def test_irsa_rejection_is_not_retried_as_an_outage(monkeypatch, status):
    import httpx

    from api.errors import ApiError
    calls = []

    def answer(method, url, **kw):
        calls.append(url)
        return httpx.Response(status, request=httpx.Request(method, url))
    monkeypatch.setattr(httpx, "request", answer)
    monkeypatch.setattr(irsa.time, "sleep", lambda s: pytest.fail("a permanent 4xx must not be retried"))
    with pytest.raises(ApiError) as e:
        irsa._send("GET", "https://irsa.example/cutout")
    assert len(calls) == 1 and str(status) in e.value.message and "отклонил" in e.value.message


@pytest.mark.parametrize("status", [408, 429, 500, 502, 503])
def test_irsa_transient_failure_keeps_retrying(monkeypatch, status):
    import httpx
    from api.errors import ApiError

    calls, pauses = [], []

    def answer(method, url, **kw):
        calls.append(url)
        return httpx.Response(status, request=httpx.Request(method, url))

    monkeypatch.setattr(httpx, "request", answer)
    monkeypatch.setattr(irsa.time, "sleep", pauses.append)
    with pytest.raises(ApiError) as error:
        irsa._send("GET", "https://irsa.example/cutout")
    assert len(calls) == 3 and pauses == [5, 10]
    assert error.value.status == 502 and error.value.code == "upstream"
    assert str(status) in error.value.message
