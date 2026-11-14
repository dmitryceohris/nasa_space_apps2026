"""B1 regressions through /api/view; only archive input is controlled."""
import time

import numpy as np
import pytest
from fastapi.testclient import TestClient

from api import irsa, jobs
from api.main import app
from test_api_pipeline import SOURCE, OVERFLOW, synthetic  # noqa: F401
from test_processing import CENTER, fake_exposure, offset


def view(mode="coadd"):
    with TestClient(app) as client:
        for _ in range(400):
            response = client.get("/api/view", params={
                "ra": CENTER.ra.deg, "dec": CENTER.dec.deg,
                "fov": 3, "det": 3, "mode": mode,
            })
            if response.status_code != 202:
                break
            time.sleep(.01)
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["status"] == "ready", result
    return result


def test_native_bad_pixel_is_removed_before_interpolation_through_http(synthetic, monkeypatch):
    import reproject
    original_load = jobs.load_exposure
    original_reproject = reproject.reproject_interp
    seen = []

    def load(*args, **kwargs):
        exp = original_load(*args, **kwargs)
        exp.image[20, 29] = 1e12
        return exp

    def project(input_data, *args, **kwargs):
        arr, _ = input_data
        if kwargs.get("order") == "bilinear":
            assert np.isnan(arr[20, 29])  # rejected in BOTH image and variance
            assert np.isfinite(arr[20, 20])  # SOURCE survives
            seen.append(True)
        return original_reproject(input_data, *args, **kwargs)

    monkeypatch.setattr(jobs, "load_exposure", load)
    monkeypatch.setattr(reproject, "reproject_interp", project)
    result = view()
    assert len(seen) == 8  # four real reprojections, image + variance
    assert all(f["measurement"]["status"] == "ok" for f in result["frames"])
    for f in result["frames"]:
        assert np.linalg.norm(np.array(f["measured_px"]) - 14) < 1
        assert len(f["measured_candidates"]) == 1


def test_only_used_exposures_contribute_time_and_spectral_samples(synthetic, monkeypatch):
    original = jobs.load_exposure

    def load(path, **kwargs):
        exp = original(path, **kwargs)
        # Metadata time differs from FITS time; u1 is downloaded but unusable.
        exp.mjd += 5
        exp.time_source = "MJD-AVG"
        if path == "u1":
            exp.flags[:] = OVERFLOW
        return exp

    monkeypatch.setattr(jobs, "load_exposure", load)
    result = view()
    first, second = result["frames"]
    assert result["exposure_counts"] == {"selected": 4, "failed": 0, "used": 3, "rejected": 1}
    assert first["mjd"] == 105 and first["mjd_span"] == [105, 105]
    assert first["obs_ids"] == ["T100.0"]
    assert first["time_sources"] == [{"obs_id": "T100.0", "source": "MJD-AVG", "assumptions": []}]
    assert [s["obs_id"] for s in first["wavelength_samples"]] == ["T100.0"]
    assert second["mjd"] == 305.5 and second["mjd_span"] == [305, 306]


def test_diff_keeps_signed_lobes_and_uses_only_positive_position(synthetic, monkeypatch):
    old, new = offset(CENTER, 0, -25), offset(CENTER, 0, 25)
    monkeypatch.setattr(jobs, "load_exposure", lambda path, **kw:
                        fake_exposure(old if path in ("u0", "u1") else new,
                                      100 if path in ("u0", "u1") else 300, seed=int(path[1:])))
    result = view("diff")
    n = result["grid"]["height_px"]
    expected_y = [(n - 1) / 2 + 25 / 6.15, (n - 1) / 2 - 25 / 6.15]
    for frame, y in zip(result["frames"], expected_y):
        positive = [c for c in frame["measured_candidates"] if c["sign"] == 1]
        negative = [c for c in frame["measured_candidates"] if c["sign"] == -1]
        assert len(positive) == len(negative) == 1
        assert positive[0]["xy"][1] == pytest.approx(y, abs=.5)
        assert negative[0]["xy"][1] == pytest.approx(n - 1 - y, abs=.5)
        if frame["measured_px"] is not None:
            assert frame["measured_px"] == positive[0]["xy"]
    assert result["processing"]["difference"]["units"] == "robust_residual_scale"


@pytest.mark.parametrize("mode", ["single", "coadd", "diff"])
def test_nan_edge_and_background_do_not_create_a_source(synthetic, monkeypatch, mode):
    def load(path, **kwargs):
        exp = fake_exposure(None, 100 if path in ("u0", "u1") else 300, seed=int(path[1:]))
        exp.image += 100 if exp.mjd > 200 else 0
        exp.image[:, :5] = np.nan
        return exp

    monkeypatch.setattr(jobs, "load_exposure", load)
    result = view(mode)
    assert all(f["measured_px"] is None and not f["measured_candidates"] for f in result["frames"])


def test_small_band_change_also_warns_against_variability(synthetic, monkeypatch):
    original = jobs.load_exposure

    def load(*args, **kwargs):
        exp = original(*args, **kwargs)
        # A change <5% can still be a different narrow LVF band.
        exp.wavelength_at_sky = lambda coord: (1.70 if exp.mjd < 200 else 1.72, .04)
        return exp

    monkeypatch.setattr(jobs, "load_exposure", load)
    result = view("diff")
    assert any("полос" in w and "переменн" in w for w in result["warnings"])


def test_api_exposes_unverified_registration_and_candidate_scope(synthetic):
    result = view()
    assert result["processing"]["registration"] == {
        "name": "header_wcs_reprojection", "independently_verified": False,
        "astrometric_precision_arcsec": None,
    }
    assert any("регистрац" in w for w in result["warnings"])
    assert any("кандидат" in w for w in result["warnings"])
    assert all(f["quality"]["verified"] is False for f in result["frames"])
