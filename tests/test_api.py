"""Offline tests: no IRSA/SIMBAD traffic. Run: .venv/bin/python -m pytest -q tests"""
import os

os.environ["PREWARM"] = "0"  # before importing the app

import pytest
from fastapi.testclient import TestClient

from api import irsa, render, resolve
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


def test_select_gap_and_even_spacing():
    ex = [dict(uri=str(i), mjd=100 + i, lam=[1, 2]) for i in range(10)] + [dict(uri="late", mjd=200, lam=[1, 2])]
    cl = irsa.select(ex, 3)
    assert len(cl) == 2 and len(cl[0]) == 3 and cl[0][0]["uri"] == "0" and cl[0][-1]["uri"] == "9"


def test_obs_id():
    assert irsa.obs_id("ibe/data/spherex/dr1/level2/2025W33_2C/x/6/level2_2025W33_2C_0386_3D6_spx_l2b-v27.fits") == "2025W33_2C_0386_3"
    assert irsa.obs_id("https://irsa.ipac.caltech.edu/datalink/links/spherex?ID=ivo://irsa.ipac/spherex_dr1?2025W38_1A_0596_1/D6") == "2025W38_1A_0596_1"


def test_barnard_proper_motion():  # SIMBAD-style record; shift over ~25 yr must be ~ pm * dt
    o = {"ra": 269.45207696, "dec": 4.69336497, "pmra_masyr": -801.551, "pmdec_masyr": 10362.394,
         "parallax_mas": 546.9759, "ref_epoch": 2000.0}
    ra, dec = resolve.propagate(o, 60903.0)
    assert 255 < (dec - o["dec"]) * 3600 < 275  # 10.36 arcsec/yr north for ~25.6 yr
