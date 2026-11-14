"""Measurement on synthetic data: NaN and bad pixels, signed residuals, ambiguous cases, a stationary source,
background subtraction and epoch limits. No network."""
import numpy as np
import pytest

from api.centroid import METHOD, find, measure
from Exposure.exposure import DEFAULT_FLAG_BITS
from processing import coadd, difference, group_epochs, make_grid
from processing.epochs import smear_px
from processing.frames import difference_details
from processing.pipeline import frames_from_exposures
from test_processing import CENTER, PIX, fake_exposure, offset


def blob(shape, x, y, s, amp):
    yy, xx = np.indices(shape)
    return amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s**2))


def noise(shape, seed=0, sigma=1.0):
    return np.random.default_rng(seed).normal(0, sigma, shape)


def dist(xy, x, y):
    return float(np.hypot(xy[0] - x, xy[1] - y))


# ------------------------------------------------------------------ NaN and bad pixels
def test_nan_strip_next_to_star_is_not_part_of_it():
    """Old behaviour: every NaN pixel joined the object with the PEAK weight, so a masked strip next to a star
    dragged the centre into the strip by many pixels."""
    img = noise((40, 40), 1) + blob((40, 40), 12.0, 20.0, 2.0, 40)
    img[:, :8] = np.nan  # masked columns 4 px from the star: its wings touch them
    m = measure(img)
    assert m["status"] == "ok" and dist(m["xy"], 12.0, 20.0) < 0.4
    assert "touches_nodata" in m["flags"] and m["quality"] == "poor"


def test_masked_half_frame_is_never_an_object():
    img = noise((40, 40), 2)
    img[:, 20:] = np.nan  # footprint covers only half of the grid
    assert find(img) == []
    img[:, :20] += blob((40, 20), 8.0, 25.0, 1.5, 30)[:, :20]
    m = measure(img)
    assert m["status"] == "ok" and dist(m["xy"], 8.0, 25.0) < 0.4 and "touches_nodata" not in m["flags"]


def test_bad_pixel_on_wing_does_not_pull_centre():
    """A masked pixel on one side of a star used to get the peak value -> pulled the centre towards it."""
    img = noise((40, 40), 3, 0.3) + blob((40, 40), 20.0, 20.0, 2.0, 40)
    img[20, 22] = np.nan
    img[21, 22] = np.inf  # inf is bad data too
    m = measure(img)
    assert m["status"] == "ok" and dist(m["xy"], 20.0, 20.0) < 0.1
    assert "masked_pixels" in m["flags"] and m["quality"] == "fair"


def test_saturated_core_labelled_when_mask_given():
    img = noise((40, 40), 4) + blob((40, 40), 22.0, 18.0, 3.0, 80)
    img[17:20, 21:24] = np.nan
    sat = np.zeros((40, 40), bool)
    sat[17:20, 21:24] = True
    m = measure(img, saturated=sat)
    assert m["status"] == "ok" and dist(m["xy"], 22.0, 18.0) < 0.3 and "saturated_core" in m["flags"]


def test_lone_hot_pixel_rejected_but_undersampled_star_kept():
    img = noise((40, 40), 5)
    img[10, 30] = 200.0  # cosmic ray / hot pixel: neighbours are pure noise
    m = measure(img)
    assert m["status"] == "none" and m["n_rejected_hot"] == 1
    img = noise((40, 40), 5) + blob((40, 40), 20.3, 19.8, 0.6, 60)  # real but barely resolved star
    m = measure(img)
    assert m["status"] == "ok" and m["n_rejected_hot"] == 0 and dist(m["xy"], 20.3, 19.8) < 0.3


def test_object_at_frame_edge_is_flagged():
    img = noise((40, 40), 6) + blob((40, 40), 1.0, 20.0, 2.0, 40)
    c = find(img)
    assert len(c) == 1 and "touches_edge" in c[0]["flags"] and c[0]["quality"] == "poor"


# ------------------------------------------------------------------ difference frames: + and - separately
def test_dipole_measured_on_positive_residual_only():
    """Barnard moves ~1 px between epochs: "+" and "-" touch. With |value| they merged into one blob centred
    between the two positions."""
    d = noise((40, 40), 7) + blob((40, 40), 17.5, 20.0, 1.5, 30) - blob((40, 40), 21.5, 20.0, 1.5, 30)
    m = measure(d, signed=True)
    assert m["status"] == "ok" and dist(m["xy"], 17.5, 20.0) < 0.6
    neg = [c for c in m["candidates"] if c["sign"] == -1]
    assert len(neg) == 1 and dist(neg[0]["xy"], 21.5, 20.0) < 0.6
    merged = measure(np.abs(d))  # what |value| gave: one blob between the two positions
    assert merged["status"] == "ok" and abs(merged["xy"][0] - 19.5) < 0.5


def test_only_negative_residual_gives_no_position():
    d = noise((40, 40), 8) - blob((40, 40), 20, 20, 2.0, 40)
    m = measure(d, signed=True)
    assert m["status"] == "none" and m["xy"] is None and "отрицательный" in m["reason"]
    assert [c["sign"] for c in m["candidates"]] == [-1]


# ------------------------------------------------------------------ ambiguous cases
def test_two_equal_stars_symmetric_about_centre_are_ambiguous():
    img = noise((41, 41), 9) + blob((41, 41), 12, 20, 1.5, 40) + blob((41, 41), 28, 20, 1.5, 40)
    m = measure(img)
    assert m["status"] == "ambiguous" and m["xy"] is None and len(m["candidates"]) >= 1


def test_largest_and_central_differ_is_ambiguous():
    img = noise((60, 60), 10) + blob((60, 60), 12, 15, 4.0, 40) + blob((60, 60), 29.5, 29.5, 1.2, 60)
    m = measure(img)
    assert m["status"] == "ambiguous" and m["xy"] is None
    assert {c["kind"] for c in m["candidates"]} == {"largest", "central"}


def test_blended_pair_is_one_elongated_object():
    img = noise((40, 40), 11) + blob((40, 40), 18, 20, 1.5, 40) + blob((40, 40), 22.5, 20, 1.5, 40)
    m = measure(img)
    assert m["status"] == "ok" and "elongated" in m["flags"] and m["quality"] != "good"
    single = measure(noise((40, 40), 11) + blob((40, 40), 20, 20, 1.5, 40))
    assert "elongated" not in single["flags"] and single["quality"] == "good"


def test_noise_blobs_never_give_a_position():
    """A 3-sigma threshold on a smoothed frame finds noise blobs in ~1 frame of 3; all of them have S/N < 5."""
    results = [measure(noise((40, 40), 1000 + s)) for s in range(100)]
    assert all(m["status"] == "none" for m in results)
    assert sum(m["n_low_snr"] for m in results) > 0  # they were seen, counted and refused
    m = next(m for m in results if m["n_low_snr"])
    assert "S/N" in m["reason"]


def test_faint_star_next_to_noise_blob_is_not_ambiguous():
    img = noise((40, 40), 12) + blob((40, 40), 20, 20, 1.5, 10)
    m = measure(img)
    assert m["status"] == "ok" and dist(m["xy"], 20, 20) < 0.5


# ------------------------------------------------------------------ error estimate
def test_error_estimate_matches_scatter():
    xs, errs = [], []
    for seed in range(40):
        m = measure(noise((40, 40), 100 + seed) + blob((40, 40), 20.2, 19.7, 1.5, 12))
        assert m["status"] == "ok"
        xs.append(m["xy"][0])
        errs.append(m["err_px"][0])
    ratio = np.std(xs) / np.mean(errs)
    assert 0.5 < ratio < 2.0, ratio
    bright = measure(noise((40, 40), 1) + blob((40, 40), 20.2, 19.7, 1.5, 120))
    assert bright["err_px"][0] < np.mean(errs) / 4


def test_method_is_described():
    assert METHOD["params"]["k_sigma"] == 3.0 and METHOD["limitations"] and METHOD["not_implemented"]
    assert measure(noise((40, 40)))["method"] == METHOD["name"]


# ------------------------------------------------------------------ artificial stationary source, end to end
DITHER = [(0.3, -0.2), (-0.35, 0.15), (0.1, 0.4)]


def _epoch(mjd0, seed, star, shifts=DITHER, flux=3.0):
    return [fake_exposure(star, mjd0 + i, flux_mjysr=flux, sigma_px=1.5, noise=0.05, seed=seed + i, shift_px=s)
            for i, s in enumerate(shifts)]


def test_stationary_source_same_centre_every_epoch_and_no_diff_signal():
    star = offset(CENTER, 2.0, -1.5)  # off the pixel grid
    wcs, shape = make_grid(CENTER, 3.0)
    tx, ty = (float(v) for v in wcs.world_to_pixel(star))
    shifts = [DITHER, DITHER[::-1], [(0.0, 0.0), (0.45, 0.45), (-0.2, -0.3)]]  # different dithers per epoch
    frames = [coadd(_epoch(60900 + 180 * k, 10 * k, star, shifts[k]), wcs, shape, CENTER) for k in range(3)]
    centres = []
    for f in frames:
        m = measure(f.image, saturated=f.contam)
        assert m["status"] == "ok" and m["quality"] == "good"
        centres.append(m["xy"])
        assert dist(m["xy"], tx, ty) < 0.15
    assert np.ptp(np.array(centres), axis=0).max() < 0.15  # no apparent motion

    dd = difference_details(frames)
    for d in dd["diffs"]:
        v = d[np.isfinite(d)]
        assert 0.8 < 1.4826 * np.median(np.abs(v - np.median(v))) < 1.2  # 1 unit = 1 sigma
        yy, xx = np.indices(shape)
        near = np.hypot(xx - tx, yy - ty) < 4
        assert np.nanmax(np.abs(d[near])) < 6  # no significant residual at the star
        assert measure(d, signed=True)["status"] != "ok" or all(
            dist(c["xy"], tx, ty) > 3 for c in measure(d, signed=True)["candidates"])


def test_moving_source_two_epochs_gives_full_amplitude_pair():
    """Leave-one-out template: with two epochs the template is the OTHER epoch, not the mean of both
    (which halved the residual and correlated the noise)."""
    old, new = offset(CENTER, 0, -15), offset(CENTER, 0, 15)
    wcs, shape = make_grid(CENTER, 3.0)
    f1 = coadd(_epoch(60900, 0, old), wcs, shape, CENTER)
    f2 = coadd(_epoch(61100, 10, new), wcs, shape, CENTER)
    (d1, d2), _ = difference([f1, f2])
    xo, yo = (float(v) for v in wcs.world_to_pixel(old))
    xn, yn = (float(v) for v in wcs.world_to_pixel(new))
    m1, m2 = measure(d1, signed=True), measure(d2, signed=True)
    # the "-" of the other epoch, 5 px away, eats into the facing wing of the "+": ~0.35 px bias, documented
    assert m1["status"] == "ok" and dist(m1["xy"], xo, yo) < 0.5
    assert m2["status"] == "ok" and dist(m2["xy"], xn, yn) < 0.5
    # antisymmetric: the "-" of one epoch is the "+" of the other
    assert np.nanmax(d1) == pytest.approx(-np.nanmin(d2), rel=0.1)


def test_diff_noise_is_unit_on_pure_noise():
    wcs, shape = make_grid(CENTER, 3.0)
    frames = [coadd(_epoch(60900 + 180 * k, 30 * k, None), wcs, shape, CENTER) for k in range(4)]
    raw = difference_details(frames, renorm=False)
    dd = difference_details(frames)
    for d, s in zip(dd["diffs"], dd["noise_scale"]):
        v = d[np.isfinite(d)]
        assert abs(np.std(v) - 1) < 0.1 and s is not None
    # bilinear reprojection smooths the noise: the propagated sigma is too large, renorm corrects it
    assert all(s < 0.95 for s in dd["noise_scale"])
    assert np.nanstd(raw["diffs"][0]) < 0.95


# ------------------------------------------------------------------ background subtraction
def test_background_offset_removed_and_source_flux_kept():
    """Zodi model misses +0.3 MJy/sr; a bright extended source flagged SOURCE covers part of the field."""
    wcs, shape = make_grid(CENTER, 6.0)
    star = offset(CENTER, 10, 10)
    exps = []
    for i in range(3):
        e = fake_exposure(star, 60900 + i, flux_mjysr=5.0, sigma_px=3.0, noise=0.05, seed=i, n=81)
        e.image += 0.3
        x0, y0 = e.sky_to_pixel(star)
        yy, xx = np.indices(e.shape)
        e.flags[np.hypot(xx - x0, yy - y0) < 10] |= 1 << DEFAULT_FLAG_BITS["SOURCE"]
        exps.append(e)
    f = coadd(exps, wcs, shape, CENTER)
    assert all(o == pytest.approx(0.3, abs=0.01) for o in f.bg_offsets)
    xs, ys = (float(v) for v in wcs.world_to_pixel(star))
    yy, xx = np.indices(shape)
    far = np.hypot(xx - xs, yy - ys) > 15
    assert abs(np.nanmedian(f.image[far])) < 0.01
    assert np.nanmax(f.image) == pytest.approx(5.0, rel=0.05)  # source not eaten by the background


def test_background_without_source_flag_still_close():
    wcs, shape = make_grid(CENTER, 6.0)
    exps = [fake_exposure(offset(CENTER, 10, 10), 60900 + i, flux_mjysr=5.0, sigma_px=3.0, noise=0.05, seed=i,
                          n=81) for i in range(2)]
    for e in exps:
        e.image += 0.3
    f = coadd(exps, wcs, shape, CENTER)
    assert all(o == pytest.approx(0.3, abs=0.01) for o in f.bg_offsets)  # sigma clipping protects it


def test_background_skipped_when_too_few_pixels():
    from processing.frames import MIN_BG_PIX, _residual_background

    img = np.full((10, 10), 0.3)
    good = np.zeros((10, 10), bool)
    good.flat[: MIN_BG_PIX - 1] = True
    assert _residual_background(img, good, np.zeros_like(good)) is None  # not enough sky: no guess
    good.flat[: MIN_BG_PIX] = True
    assert _residual_background(img, good, good.copy()) == pytest.approx(0.3)  # all SOURCE -> falls back


# ------------------------------------------------------------------ epoch limits
def test_epoch_limits_hold_on_random_cadence():
    rng = np.random.default_rng(0)
    for _ in range(50):
        mjds = np.sort(np.concatenate([rng.uniform(0, 90, rng.integers(1, 40)),
                                       rng.uniform(180, 200, rng.integers(0, 10))]))
        rng.shuffle(mjds)
        groups = group_epochs(mjds, gap_days=30, max_span_days=20)
        flat = sorted(i for g in groups for i in g)
        assert flat == list(range(len(mjds)))  # every exposure exactly once
        for g in groups:
            assert np.ptp(mjds[g]) <= 20
        starts = [mjds[g].min() for g in groups]
        assert starts == sorted(starts)
        for a, b in zip(groups[:-1], groups[1:]):
            assert mjds[a].max() <= mjds[b].min()  # epochs do not overlap in time


def test_epoch_limits_bad_parameters():
    with pytest.raises(ValueError):
        group_epochs([1, 2], gap_days=0)
    with pytest.raises(ValueError):
        group_epochs([1, 2], max_span_days=-1)
    with pytest.raises(ValueError):
        group_epochs([1, np.nan])
    assert group_epochs([5, 5, 5], max_span_days=0) == [[0, 1, 2]]  # same instant is one epoch


def test_continuous_coverage_frames_respect_span():
    exps = [fake_exposure(CENTER, 60900 + 2 * i, seed=i) for i in range(30)]  # 58 days, ecliptic-pole style
    frames, _, _ = frames_from_exposures(exps, CENTER, 2.0)
    assert len(frames) == 3 and all(f.span_days <= 20 for f in frames)
    assert sum(f.n_exp for f in frames) == 30


def test_smear_within_epoch():
    # Barnard: 10.4"/yr, 20-day epoch -> 0.09 px; an asteroid-like 1000"/yr -> ~9 px
    assert smear_px(10.4, 20, PIX) < 0.1
    assert smear_px(1000, 20, PIX) > 5
