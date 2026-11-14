"""Тесты API без сети: тяжёлые функции подменяются заглушками."""
import importlib
import os

import pytest


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("PREWARM", "0")
    import api.jobs
    import api.main
    importlib.reload(api.jobs)
    main = importlib.reload(api.main)

    calls = {"view": 0, "spectrum": 0}

    def fake_view(center, fov, det, n, out_dir, prefix, target=None, progress=None):
        calls["view"] += 1
        progress("download", 1, 2, "Скачано 1 из 2")
        (out_dir / "f0.png").write_bytes(b"png")
        return {"status": "ready", "frames": [{"index": 0, "urls": {"coadd": f"{prefix}/f0.png"},
                                               "predicted_px": [1, 2]}], "warnings": []}

    def fake_spectrum(center, coord_at, per_epoch=2, progress=None):
        calls["spectrum"] += 1
        coord_at(61000.0)
        pts = [{"wave_um": w, "flux_ujy": 100.0 + i, "err_ujy": 1.0, "mjd": 61000.0 + i, "det": 2,
                "ok": True, "flags": []} for i, w in enumerate([1.30, 1.35, 1.40, 2.0])]
        return {"status": "ready", "points": pts, "n_points": len(pts), "n_saturated": 0}

    monkeypatch.setattr(main, "build_view", fake_view)
    monkeypatch.setattr(main, "build_spectrum", fake_spectrum)
    from fastapi.testclient import TestClient
    with TestClient(main.app) as c:
        c.calls, c.jobs = calls, main.jobs
        yield c


def test_health_and_targets(client):
    assert client.get("/api/health").json()["status"] == "ok"
    ids = {t["id"] for t in client.get("/api/targets").json()}
    assert ids == {"barnard", "epsindib", "wise0855", "rwaur", "rdra"}


def test_bad_params(client):
    r = client.get("/api/view", params={"ra": 10, "dec": 100})
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"
    assert client.get("/api/view").status_code == 400
    assert client.get("/api/view", params={"ra": 10}).status_code == 400
    assert client.get("/api/view", params={"target": "barnard", "det": 7}).status_code == 400
    assert client.get("/api/view", params={"target": "barnard", "fov": 50}).status_code == 400


def test_resolve_preset(client):
    r = client.get("/api/resolve", params={"target": "Barnard's Star"}).json()
    assert r["id"] == "barnard" and r["dec_now"] > r["dec"]   # звезда уходит на север


def test_view_202_then_200_and_cache(client):
    p = {"target": "wise0855", "det": 6, "fov": 3}
    r1 = client.get("/api/view", params=p)
    assert r1.status_code in (200, 202)
    client.jobs.wait_all()
    r2 = client.get("/api/view", params=p)
    assert r2.status_code == 200 and r2.json()["status"] == "ready"
    q = r2.json()["query"]
    assert q["det"] == 6 and q["target"]["id"] == "wise0855" and q["band_um"] == [4.42, 5.0]
    png = r2.json()["frames"][0]["urls"]["coadd"]
    assert client.get(png).status_code == 200                  # картинка отдаётся статикой
    assert client.get("/api/view", params=p).status_code == 200
    assert client.calls["view"] == 1                           # второй раз — из кэша


def test_sexagesimal_coordinates(client):
    r = client.get("/api/view", params={"ra": "08h55m10.831s", "dec": "-7d14m42.53s"})
    assert r.status_code in (200, 202)
    client.jobs.wait_all()
    q = client.get("/api/view", params={"ra": "08h55m10.831s", "dec": "-7d14m42.53s"}).json()["query"]
    assert q["ra"] == pytest.approx(133.79513, abs=1e-4) and q["dec"] == pytest.approx(-7.24515, abs=1e-4)


def test_spectrum_and_lightcurve(client):
    p = {"target": "rwaur"}
    client.get("/api/spectrum", params=p)
    client.jobs.wait_all()
    s = client.get("/api/spectrum", params=p).json()
    assert s["n_points"] == 4
    lc = client.get("/api/lightcurve", params={**p, "det": 2}).json()
    assert lc["status"] == "ready"
    assert [round(x["wave_um"], 2) for x in lc["points"]] == [1.30, 1.35, 1.40]   # 2.0 мкм вне окна
    assert client.calls["spectrum"] == 1


def test_frontend_served(client):
    r = client.get("/")
    assert r.status_code == 200 and "SPHEREx Sky Viewer" in r.text
