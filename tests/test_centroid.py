"""Synthetic frames: noise, different object sizes, saturation hole. Run with PYTHONPATH set to the project root."""
import numpy as np

from api.centroid import find


def blob(shape, x, y, sx, amp):
    yy, xx = np.indices(shape)
    return amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sx**2))


def noisy(shape, seed=0, sigma=1.0):
    return np.random.default_rng(seed).normal(0, sigma, shape)


def close(xy, x, y, tol=0.7):
    return abs(xy[0] - x) < tol and abs(xy[1] - y) < tol


def test_single_object_one_centre():
    img = noisy((40, 40)) + blob((40, 40), 20.3, 19.6, 2.0, 40)
    c = find(img)
    assert len(c) == 1 and c[0]["kind"] == "largest+central" and close(c[0]["xy"], 20.3, 19.6)


def test_large_offcentre_and_small_central_give_two_centres():
    img = noisy((60, 60), 1) + blob((60, 60), 12, 15, 4.0, 40) + blob((60, 60), 29.5, 29.5, 1.2, 60)
    c = {x["kind"]: x for x in find(img)}
    assert set(c) == {"largest", "central"}
    assert close(c["largest"]["xy"], 12, 15) and close(c["central"]["xy"], 29.5, 29.5)


def test_other_scale_same_answer():
    for sx in (1.0, 2.5, 5.0):
        img = noisy((80, 80), 2) + blob((80, 80), 41.2, 38.4, sx, 60)
        c = find(img)
        assert len(c) == 1 and close(c[0]["xy"], 41.2, 38.4), sx


def test_saturation_hole_does_not_bias():
    img = noisy((40, 40), 3) + blob((40, 40), 22, 18, 3.0, 80)
    img[17:20, 21:24] = np.nan  # masked core
    c = find(img)
    assert len(c) == 1 and close(c[0]["xy"], 22, 18, 0.5)


def test_pure_noise_and_flat_give_nothing():
    assert find(noisy((40, 40), 4)) == []
    assert find(np.ones((40, 40))) == []
    assert find(np.full((40, 40), np.nan)) == []


def test_signed_diff_counts_dark_residual():
    img = noisy((40, 40), 5) - blob((40, 40), 15, 25, 2.0, 50)
    assert len(find(img, signed=True)) == 1 and find(img) == []
