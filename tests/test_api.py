"""Offline tests: no IRSA/SIMBAD traffic. Run: .venv/bin/python -m pytest -q tests"""
import os
import time

os.environ["PREWARM"] = "0"  # before importing the app

import pytest
from fastapi.testclient import TestClient

from Exposure import DETECTOR_BANDS
from api import cache, irsa, jobs, render, resolve
from api.main import app
from api.params import parse

c = TestClient(app)


def test_health():
    assert c.get("/api/health").json() == {"status": "ok", "version": "0.1.0"}


def test_targets_are_the_spec_presets():
    t = c.get("/api/targets").json()
    assert [x["id"] for x in t] == ["barnard", "epsindib", "wise0855", "rwaur", "rdra"]


@pytest.mark.parametrize("q", ["dec=100&ra=10", "ra=10", "", "target=x&fov=50", "target=x&det=9",
                               "target=x&mode=zzz", "target=x&max_per_epoch=0"])
def test_bad_param(q):
    r = c.get("/api/view?" + q)
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"


def test_dec_message():
    assert c.get("/api/view?ra=10&dec=100").json()["error"]["message"] == "dec должен быть в диапазоне -90..90"


def test_unknown_name_is_404(monkeypatch):
    monkeypatch.setattr(resolve, "_fetch_sesame", lambda name: ["#=Sc=Simbad\n#! nothing\n"])
    r = c.get("/api/resolve?target=asdfgh")
    assert r.status_code == 404 and r.json()["error"]["code"] == "not_found"


def test_all_presets_resolve_from_team_table_without_network(monkeypatch):
    def no_network(name):
        raise AssertionError("SIMBAD must not be asked for the five presets")

    monkeypatch.setattr(resolve, "_fetch_sesame", no_network)
    for tid in ["barnard", "epsindib", "wise0855", "rwaur", "rdra"]:
        o = resolve.resolve(tid)
        assert o["id"] == tid and o["pmra_masyr"] is not None and o["source"]
    assert resolve.resolve("WISE 0855-0714")["id"] == "wise0855"  # the spec name that SIMBAD does not know
    assert c.get("/api/resolve?target=wise0855").status_code == 200


def test_spectrum_needs_target_or_coords():
    r = c.get("/api/spectrum")
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"


@pytest.mark.parametrize("q", ["wave_min=3&wave_max=2", "wave_min=abc", "wave_min=0&wave_max=1"])
def test_lightcurve_bad_window(q):
    r = c.get("/api/lightcurve?ra=10&dec=10&" + q)
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"


def test_default_window_is_5_percent_around_band_middle():
    lo, hi, _ = DETECTOR_BANDS[3]
    wlo, whi = jobs.default_window(3)
    assert wlo == pytest.approx((lo + hi) / 2 * 0.95) and whi == pytest.approx((lo + hi) / 2 * 1.05)


def test_spectrum_job_202_then_200_and_lightcurve(monkeypatch, tmp_path):
    """build_spectrum is replaced by a fake: only the job mechanism and the light-curve window are tested."""
    monkeypatch.setattr(cache, "CACHE", tmp_path)

    def pt(w, mjd):
        return {"wave_um": w, "bw_um": 0.02, "flux_ujy": 5.0, "err_ujy": 1.0, "mjd": mjd, "det": 3,
                "obs_id": "a", "ok": True, "flags": []}

    def fake_build(center, coord_at, per_epoch=2, progress=None):
        progress("download", 1, 2, "fake")
        time.sleep(0.2)
        return {"status": "ready", "points": [pt(2.0, 60000.0), pt(2.1, 60100.0), pt(4.9, 60200.0)],
                "n_points": 3, "n_saturated": 0, "credits": "c", "bands": {}}

    monkeypatch.setattr(jobs, "build_spectrum", fake_build)
    first = c.get("/api/spectrum?ra=33.3&dec=44.4")
    assert first.status_code == 202 and first.json()["status"] == "running"
    for _ in range(100):
        r = c.get("/api/spectrum?ra=33.3&dec=44.4")
        if r.status_code == 200:
            break
        time.sleep(0.05)
    assert r.status_code == 200 and r.json()["n_points"] == 3 and "v" not in r.json()
    lc = c.get("/api/lightcurve?ra=33.3&dec=44.4&wave_min=1.9&wave_max=2.2").json()
    assert lc["wave_window_um"] == [1.9, 2.2] and [p["mjd"] for p in lc["points"]] == [60000.0, 60100.0]
    # flat arrays the frontend plots directly (frontend/index.html: d.dates / d.flux / d.flux_err)
    assert lc["dates"] == ["2023-02-25", "2023-06-05"] and lc["flux"] == [5.0, 5.0] and lc["flux_err"] == [1.0, 1.0]


def test_spectrum_flat_arrays_keep_only_reliable_points():
    from processing.pipeline import _good

    pts = [{"wave_um": 1.0, "flux_ujy": 2.0, "err_ujy": 0.1, "ok": True},
           {"wave_um": 1.1, "flux_ujy": 9.0, "err_ujy": 0.1, "ok": False},   # saturated: not plotted
           {"wave_um": 1.2, "flux_ujy": None, "err_ujy": None, "ok": True}]
    assert [p["wave_um"] for p in _good(pts)] == [1.0]


def test_root_serves_team_frontend_and_dev_page():
    r = c.get("/")
    assert r.status_code == 200 and "SPHEREx Sky Viewer" in r.text and "pollJson" in r.text
    assert c.get("/dev/").status_code == 200


def test_export_before_ready_is_409():
    assert c.get("/api/export.gif?ra=11.1&dec=22.2&fov=7").status_code == 409


def test_coords_degrees_and_sexagesimal():
    a = parse(ra="133.795", dec="-7.245")
    b = parse(ra="08h55m10.8s", dec="-07d14m42s")
    assert a.ra == 133.795 and abs(b.ra - 133.795) < 0.01 and abs(b.dec + 7.245) < 0.01


def test_radec_beats_target():
    p = parse(target="barnard", ra="10", dec="20")
    assert (p.ra, p.dec) == (10.0, 20.0)


def test_png_y_flip():
    assert render.png_xy(3, 0, 10) == [3.0, 9.0]  # FITS bottom row -> PNG bottom row


def test_obs_id():
    assert irsa.obs_id("ibe/data/spherex/dr1/level2/2025W33_2C/x/6/level2_2025W33_2C_0386_3D6_spx_l2b-v27.fits") == "2025W33_2C_0386_3"
    assert irsa.obs_id("https://irsa.ipac.caltech.edu/datalink/links/spherex?ID=ivo://irsa.ipac/spherex_dr1?2025W38_1A_0596_1/D6") == "2025W38_1A_0596_1"


def test_barnard_proper_motion():  # SIMBAD-style record; shift over ~25 yr must be ~ pm * dt
    o = {"ra": 269.45207696, "dec": 4.69336497, "pmra_masyr": -801.551, "pmdec_masyr": 10362.394,
         "parallax_mas": 546.9759, "ref_epoch": 2000.0}
    ra, dec = resolve.propagate(o, 60903.0)
    assert 255 < (dec - o["dec"]) * 3600 < 275  # 10.36 arcsec/yr north for ~25.6 yr


@pytest.mark.parametrize(('ra', 'dec'), [
    ('08h55m10.831s', '-7d14m42.53s'),
    ('08:55:10.831', '-07:14:42.53'),
    ('133.79512916666667', '-7d14m42.53s'),
    ('08:55:10.831', '-7.245147222222222'),
])
def test_wise_coordinate_formats(ra, dec):
    p = parse(ra=f'  {ra} ', dec=f' {dec}  ')
    assert p.ra == pytest.approx(133.79512916666667)
    assert p.dec == pytest.approx(-7.245147222222222)


@pytest.mark.parametrize(('ra', 'dec', 'expected_ra'), [('0', '-90', 0), ('360', '90', 0), ('24:00:00', '0', 0)])
def test_coordinate_boundaries(ra, dec, expected_ra):
    p = parse(ra=ra, dec=dec)
    assert p.ra == expected_ra
    assert p.dec == float(dec)


@pytest.mark.parametrize(('ra', 'dec'), [
    ('361', '0'), ('-1', '0'), ('24:00:01', '0'), ('0', '90.01'),
    ('0', '-90:00:01'), ('nan', '0'), ('inf', '0'), ('0', '-inf'),
    ('bad', '0'), ('1:60:00', '0'), ('0', ' '), (' ', '0'),
])
def test_invalid_coordinates_do_not_start_jobs(monkeypatch, ra, dec):
    def fail(*args, **kwargs):
        pytest.fail('invalid coordinates started a job')
    monkeypatch.setattr('api.jobs.view', fail)
    r = c.get('/api/view', params={'ra': ra, 'dec': dec, 'target': 'barnard'})
    assert r.status_code == 400 and r.json()['error']['code'] == 'bad_param'


def test_explicit_coordinates_never_resolve_or_propagate_target(monkeypatch):
    def fail(*args, **kwargs):
        pytest.fail('explicit coordinates should be fixed')
    monkeypatch.setattr(resolve, 'resolve', fail)
    monkeypatch.setattr(resolve, 'propagate', fail)
    result = resolve.for_params(parse(target='barnard', ra='360', dec='-7:14:42.53'))
    assert result['center'] == pytest.approx((0, -7.245147222222222))
    assert result['motion'] is None


def test_equivalent_coordinate_requests_reach_same_job(monkeypatch):
    seen = []
    def capture(p):
        seen.append(p)
        return 200, {'status': 'ready'}
    monkeypatch.setattr('api.jobs.view', capture)
    for ra, dec in [('08:55:10.831', '-7d14m42.53s'), ('133.79512916666667', '-7.245147222222222')]:
        assert c.get('/api/view', params={'ra': ra, 'dec': dec}).status_code == 200
    assert seen[0].ra == pytest.approx(seen[1].ra)
    assert seen[0].dec == pytest.approx(seen[1].dec)
    assert (seen[0].fov, seen[0].det, seen[0].mode) == (seen[1].fov, seen[1].det, seen[1].mode)


def test_framework_errors_use_the_same_error_shape():
    r = c.get("/api/resolve")  # missing required query param
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"
    r = c.get("/api/no-such-endpoint")
    assert r.status_code == 404 and r.json()["error"]["code"] == "not_found"
