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


def test_data_states_distinguish_empty_failed_and_unusable(synthetic, monkeypatch):
    monkeypatch.setattr(irsa, 'exposures', lambda *args: [])
    empty, _ = run('coadd')
    assert empty['data_status'] == 'no_archive_results'
    assert 'никогда' not in ' '.join(empty['warnings'])
    monkeypatch.setattr(irsa, 'exposures', lambda *args: [{'uri': 'broken', 'mjd': 60000}])
    def fail(*args):
        raise OSError('controlled download failure')
    monkeypatch.setattr(irsa, 'cutout', fail)
    failed, _ = run('coadd')
    assert failed['data_status'] == 'download_failed'
    monkeypatch.setattr(irsa, 'cutout', lambda *args: 'broken')
    unusable, _ = run('coadd')
    assert unusable['data_status'] == 'unusable_data'


def test_available_data_reports_counts(synthetic):
    result, _ = run('single')
    assert result['data_status'] == 'available'
    assert result['exposure_counts'] == {'selected': 2, 'failed': 0, 'used': 2, 'rejected': 0}
    assert all(f['quality']['coverage_fraction'] > 0 for f in result['frames'])
    assert all(f['quality']['verified'] is False for f in result['frames'])


def test_poor_centroid_is_only_a_candidate(synthetic, monkeypatch):
    original = jobs.centroid.measure
    def poor(*args, **kwargs):
        result = original(*args, **kwargs)
        result['quality'] = 'poor'
        return result
    monkeypatch.setattr(jobs.centroid, 'measure', poor)
    result, _ = run('single')
    assert all(frame['measured_px'] is None for frame in result['frames'])
    assert all(frame['measured_candidates'] for frame in result['frames'])


@pytest.mark.parametrize('mode', ['single', 'coadd', 'diff', 'trail'])
def test_real_http_path_uses_processing_and_preserves_quality(synthetic, mode, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from api.main import app
    mount = next(route.app for route in app.routes if getattr(route, 'name', '') == 'cache')
    monkeypatch.setattr(mount, 'all_directories', [str(synthetic)])
    with TestClient(app) as client:
        query = {'ra': CENTER.ra.deg, 'dec': CENTER.dec.deg, 'fov': 3, 'det': 3, 'mode': mode}
        for _ in range(200):
            response = client.get('/api/view', params=query)
            if response.status_code != 202:
                break
            time.sleep(0.01)
        assert response.status_code == 200
        result = response.json()
        assert result['status'] == 'ready' and len(result['frames']) == 2
        assert result['target_info'] is None
        assert result['frames'][0]['quality']['verified'] is False
        assert result['frames'][0]['wavelength_samples'][0]['source'] == 'approximate_detector_position'
        image = client.get(result['frames'][0]['image_url'])
        assert image.status_code == 200 and image.content.startswith(b'\x89PNG')


def _http_result(mode):
    import time
    from fastapi.testclient import TestClient
    from api.main import app
    client = TestClient(app)
    for _ in range(200):
        response = client.get('/api/view', params={'ra':CENTER.ra.deg,'dec':CENTER.dec.deg,'fov':3,'det':3,'mode':mode})
        if response.status_code != 202:
            break
        time.sleep(.01)
    assert response.status_code == 200
    return response.json()


@pytest.mark.parametrize('mode', ['single','coadd','diff','trail'])
def test_one_epoch_through_http_has_no_fake_animation(synthetic, monkeypatch, mode):
    monkeypatch.setattr(irsa, 'exposures', lambda *args: [{'uri':'u0','mjd':100}])
    result = _http_result(mode)
    assert result['status'] == 'ready' and result['data_status'] == 'single_epoch'
    assert len(result['frames']) == (0 if mode == 'diff' else 1)
    assert result['trail_url'] is None
    assert any('одна эпоха' in message for message in result['warnings'])


def test_partial_http_result_retries_failed_exposure_after_short_ttl(synthetic, monkeypatch):
    import os
    import time
    def sometimes_broken(uri, *args):
        if uri == 'u0':
            raise OSError('controlled outage')
        return uri
    monkeypatch.setattr(irsa, 'cutout', sometimes_broken)
    first = _http_result('coadd')
    assert first['data_status'] == 'partial' and first['exposure_counts']['failed'] == 1
    key = first['frames'][0]['image_url'].split('/')[2]
    stamp = time.time()-61
    os.utime(synthetic/key/'result.json', (stamp,stamp))
    monkeypatch.setattr(irsa, 'cutout', lambda uri,*args: uri)
    recovered = _http_result('coadd')
    assert recovered['data_status'] == 'available'
    assert recovered['exposure_counts']['used'] == 4 and recovered['exposure_counts']['failed'] == 0


def test_stationary_star_and_changed_background_through_http(synthetic, monkeypatch):
    original = jobs.load_exposure
    def changed(*args, **kwargs):
        exp = original(*args, **kwargs)
        if exp.mjd > 200:
            exp.image = exp.image + 40
        exp.image[:, 0] = np.nan  # no-data edge must not become a source
        return exp
    monkeypatch.setattr(jobs, 'load_exposure', changed)
    result = _http_result('diff')
    assert len(result['frames']) == 2
    assert all(f['measured_px'] is None for f in result['frames'])
    assert result['processing']['difference']['units'] == 'robust_residual_scale'


def test_outside_prediction_is_hidden_through_http(synthetic, monkeypatch):
    from api import resolve
    from astropy.coordinates import SkyCoord
    # A controlled moving-target adapter: query centre is covered, old-date prediction is outside.
    monkeypatch.setattr(resolve, 'for_params', lambda p: {
        'center':(CENTER.ra.deg,CENTER.dec.deg), 'label':'controlled moving target',
        'motion':{'name':'controlled moving target','source':'fixture','pmra_masyr':0,'pmdec_masyr':0}})
    monkeypatch.setattr(jobs, '_sky', lambda r,mjd: SkyCoord(151,20,unit='deg'))
    result = _http_result('single')
    assert len(result['frames']) == 2
    assert all(frame['predicted_px'] is None and frame['prediction_status']=='outside_field' for frame in result['frames'])
    assert any('маркер скрыт' in warning for warning in result['warnings'])
