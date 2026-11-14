"""API acceptance review (docs/reviews/api-acceptance.md): each test pins one defect found in eaca1e5."""
import json
import os

os.environ["PREWARM"] = "0"

import pytest
from fastapi.testclient import TestClient

from api import cache, jobs, resolve
from api.errors import ApiError
from api.main import app

c = TestClient(app)


def test_framework_errors_use_the_documented_error_shape():
    # docs/API.md: "Все эндпоинты возвращают ошибки в одном виде"; the frontend reads d.error.message.
    r = c.get("/api/resolve")
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"
    r = c.get("/api/no-such-endpoint")
    assert r.status_code == 404 and r.json()["error"]["code"] == "not_found"


def test_finished_failed_jobs_do_not_fill_the_queue(monkeypatch):
    # Failed jobs are retained for JOB_RETENTION_SECONDS; they must not block new work as "busy".
    monkeypatch.setattr(cache, "read", lambda key: None)
    monkeypatch.setattr(jobs, "JOBS", {"a": jobs.Job("a", "1", status="error", finished=jobs.time.monotonic()),
                                       "b": jobs.Job("b", "2", status="error", finished=jobs.time.monotonic())})
    monkeypatch.setenv("MAX_PENDING_JOBS", "2")
    monkeypatch.setattr(jobs.POOL, "submit", lambda *args: None)
    status, payload = jobs._poll("new", lambda job: None)
    assert status == 202 and payload["status"] == "running"


def test_failed_downloads_without_usable_points_are_retried_soon(monkeypatch, tmp_path):
    # Spectrum: some downloads failed and the rest gave no reliable point -> state unusable_data with points.
    monkeypatch.setattr(cache, "CACHE", tmp_path)
    monkeypatch.setattr(cache.time, "time", lambda: 1000)
    cache.write("spec", {"status": "ready", "data_status": "unusable_data", "points": [{"ok": False}],
                         "exposure_counts": {"selected": 4, "downloaded": 1}})
    os.utime(tmp_path / "spec" / "result.json", (900, 900))
    assert cache.read("spec") is None


def test_interrupted_name_cache_does_not_break_name_resolution(monkeypatch, tmp_path):
    path = tmp_path / "resolve.json"
    path.write_text('{"vega": {"ra": 1', encoding="utf-8")  # interrupted write
    monkeypatch.setattr(resolve, "_FILE", path)
    monkeypatch.setattr(resolve, "_mem", {})
    sesame = "#=Sc=Simbad\n%J 279.23473479 +38.78368896 = 18:36:56.33 +38:47:01.2\n%I.0 * alf Lyr\n"
    monkeypatch.setattr(resolve, "_fetch_sesame", lambda name: [sesame])
    assert resolve.resolve("Vega")["ra"] == pytest.approx(279.23473479)
    assert "vega" in json.loads(path.read_text(encoding="utf-8"))
    assert list(tmp_path.glob("*.tmp")) == []


def test_empty_lightcurve_window_explains_itself_and_keeps_spectrum_warnings(monkeypatch):
    spec = {"status": "ready", "warnings": ["Длина волны части точек приблизительна"], "query": {}, "target_info": None,
            "points": [{"mjd": 61000.0, "wave_um": 4.0, "bw_um": 0.04, "flux_ujy": 1.0, "err_ujy": 0.1, "ok": True}]}
    monkeypatch.setattr(jobs, "spectrum_view", lambda p: (200, spec))
    status, res = jobs.lightcurve_view(None, 1.9, 2.1)
    assert status == 200 and res["flux"] == []
    assert res["warnings"][0] == spec["warnings"][0] and len(res["warnings"]) == 2
    assert spec["warnings"] == ["Длина волны части точек приблизительна"]  # cached spectrum not mutated
