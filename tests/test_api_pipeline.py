"""Frame pipeline of the API on synthetic exposures (no network): flags layer, modes, centres, trail."""
import os

os.environ["PREWARM"] = "0"

import numpy as np
import pytest
from PIL import Image

from Exposure.exposure import DEFAULT_FLAG_BITS
from api import cache, irsa, jobs
from api.params import Params
from test_processing import CENTER, fake_exposure

SOURCE, OVERFLOW = 1 << DEFAULT_FLAG_BITS["SOURCE"], 1 << DEFAULT_FLAG_BITS["OVERFLOW"]


@pytest.fixture
def synthetic(monkeypatch, tmp_path):
    """Two epochs x two exposures of one star at the field centre; SOURCE flag on the star, OVERFLOW elsewhere."""
    monkeypatch.setattr(cache, "CACHE", tmp_path)
    expos, meta = {}, []
    for i, mjd in enumerate([100.0, 101.0, 300.0, 301.0]):
        e = fake_exposure(CENTER, mjd, seed=i)
        e.flags[19:22, 19:22] = SOURCE  # the star itself
        e.flags[19:22, 28:31] = OVERFLOW  # 8 px to the side
        expos[f"u{i}"] = e
        meta.append({"uri": f"u{i}", "mjd": mjd, "lam": [2.0, 2.1]})
    monkeypatch.setattr(irsa, "exposures", lambda ra, dec, det: meta)
    monkeypatch.setattr(irsa, "cutout", lambda uri, ra, dec, size: uri)
    monkeypatch.setattr(jobs, "load_exposure", lambda path, mjd_hint=None, obs_id_hint="": expos[path])
    return tmp_path


def run(mode):
    r = {"center": (CENTER.ra.deg, CENTER.dec.deg), "motion": None, "label": "synthetic"}
    p = Params(None, CENTER.ra.deg, CENTER.dec.deg, 3.0, 3, mode, 6)
    return jobs._pipeline(jobs.Job("k" + mode, "id"), p, r), r


def alpha(path_url, tmp):
    return np.asarray(Image.open(tmp / path_url.removeprefix("/cache/")).convert("RGBA"))[..., 3]


def test_flags_layer_marks_bad_pixels_but_never_the_source(synthetic):
    res, _ = run("coadd")
    assert len(res["frames"]) == 2 and res["trail_url"] is None
    a = alpha(res["frames"][0]["flags_url"], synthetic)
    h = a.shape[0]
    assert a[h // 2, h // 2] == 0  # the star (SOURCE) is not painted
    assert a[h // 2, h // 2 + 8] > 0  # OVERFLOW is painted


def test_one_object_gives_one_centre_at_the_target(synthetic):
    res, _ = run("coadd")
    f = res["frames"][0]
    assert [c["kind"] for c in f["measured_candidates"]] == ["largest+central"]
    x, y = f["measured_px"]
    n = res["grid"]["width_px"]
    assert abs(x - (n - 1) / 2) < 1.0 and abs(y - (n - 1) / 2) < 1.0
    assert f["predicted_px"] is None  # explicit coordinates: no proper motion, nothing predicted


def test_diff_mode_gives_trail_and_signed_frames(synthetic):
    res, _ = run("diff")
    assert len(res["frames"]) == 2 and res["trail_url"]
    assert (synthetic / res["trail_url"].removeprefix("/cache/")).exists()
    # a stationary star: no "+" residual worth a position, and the diff noise was renormalised
    for f in res["frames"]:
        assert f["measurement"]["status"] != "ok"
        assert f["diff_noise_scale"] is not None
        assert all(c["sign"] in (1, -1) for c in f["measured_candidates"])
    assert res["processing"]["difference"]["name"] == "leave_one_out_median"


def test_method_quality_and_limits_are_in_the_answer(synthetic):
    res, _ = run("coadd")
    pr = res["processing"]
    assert pr["measurement"]["name"] == "threshold_moment_centroid" and pr["measurement"]["limitations"]
    assert "PSF-фиттинг" in pr["measurement"]["not_implemented"]
    assert pr["epochs"]["max_span_days"] == 20 and pr["background"]["limitations"]
    assert pr["difference"] is None  # coadd mode: no difference frames
    f = res["frames"][0]
    m = f["measurement"]
    assert m["status"] == "ok" and m["quality"] in ("good", "fair") and max(m["err_px"]) < 0.1
    assert f["span_days"] == 1.0 and len(f["bg_offsets"]) == 2
    c = f["measured_candidates"][0]
    assert {"snr", "err_px", "quality", "flags", "sign"} <= set(c) and c["sign"] is None


def test_trail_mode_shows_coadd_frames_plus_trail(synthetic):
    res, _ = run("trail")
    assert len(res["frames"]) == 2 and res["trail_url"] and res["frames"][0]["image_url"].endswith("f0.png")


def test_single_mode_uses_one_exposure_per_epoch(synthetic):
    res, _ = run("single")
    assert [f["n_exp"] for f in res["frames"]] == [1, 1]


def test_diff_needs_two_epochs(synthetic, monkeypatch):
    meta = [{"uri": "u0", "mjd": 100.0, "lam": [2.0, 2.1]}]
    monkeypatch.setattr(irsa, "exposures", lambda ra, dec, det: meta)
    res, _ = run("diff")
    assert res["frames"] == [] and any("минимум 2 эпохи" in w for w in res["warnings"])
