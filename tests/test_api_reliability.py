"""Reliability, deduplication, bounded queue, atomic cache, TTL and stale fallback tests."""
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

os.environ["PREWARM"] = "0"

import httpx
import pytest
from fastapi.testclient import TestClient

from api import cache, irsa, jobs, resolve
from api.errors import ApiError
from api.main import app
from test_processing import CENTER, fake_exposure

c = TestClient(app)


@pytest.fixture(autouse=True)
def isolated_cache(tmp_path, monkeypatch):
    """Never touch the project cache directory during tests."""
    monkeypatch.setattr(cache, "CACHE", tmp_path)
    monkeypatch.setattr(irsa, "CACHE", tmp_path)
    monkeypatch.setattr(resolve, "_FILE", tmp_path / "resolve.json")
    resolve._mem.clear()
    with jobs.LOCK:
        jobs.JOBS.clear()
    return tmp_path


def test_concurrent_deduplication_empty_and_error(monkeypatch):
    """Parallel identical requests for empty or failing jobs must invoke the runner exactly once."""
    calls = {"view": 0}
    barrier = threading.Barrier(4)

    def empty_exposures(ra, dec, det):
        calls["view"] += 1
        time.sleep(0.15)
        return []

    monkeypatch.setattr(irsa, "exposures", empty_exposures)

    def worker():
        barrier.wait()
        for _ in range(30):
            r = c.get("/api/view?ra=12.34&dec=56.78&fov=3&det=3&mode=coadd")
            if r.status_code == 200:
                return r.json()
            time.sleep(0.02)
        raise AssertionError("timed out")

    with ThreadPoolExecutor(4) as ex:
        results = [f.result() for f in [ex.submit(worker) for _ in range(4)]]

    assert all(r["status"] == "ready" and r["frames"] == [] for r in results)
    assert calls["view"] == 1, f"Expected 1 runner invocation, got {calls['view']}"


def test_concurrent_deduplication_resolve_and_irsa_cutout(monkeypatch, isolated_cache):
    """Parallel /api/resolve and parallel irsa.cutout for the same URI must do only 1 network fetch."""
    sesame_calls = {"n": 0}

    def slow_sesame(name):
        sesame_calls["n"] += 1
        time.sleep(0.15)
        return ["#=Sc=Simbad\n%J 10.0 20.0\n%P 1.0 2.0\n%X 5.0\n"]

    monkeypatch.setattr(resolve, "_fetch_sesame", slow_sesame)
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(lambda _: c.get("/api/resolve?target=VegaTest").json(), range(4)))
    assert all(r["ra"] == 10.0 for r in res)
    assert sesame_calls["n"] == 1, f"Expected 1 Sesame call, got {sesame_calls['n']}"

    http_calls = {"n": 0}

    class DummyResp:
        status_code = 200
        content = b"SIMPLE  =                    T"
        text = "uri,mjd,lo,hi\nlevel2_2025W10_1A_0001_1D3_spx.fits,60900.0,1e-6,2e-6\n"

    def slow_send(method, url, **kw):
        http_calls["n"] += 1
        time.sleep(0.15)
        return DummyResp()

    monkeypatch.setattr(irsa, "_send", slow_send)
    with ThreadPoolExecutor(4) as ex:
        paths = list(ex.map(lambda _: irsa.cutout("ibe/level2_test.fits", 10.0, 20.0, 0.1), range(4)))
    assert len(set(paths)) == 1
    assert http_calls["n"] == 1, f"Expected 1 cutout download, got {http_calls['n']}"


def test_bounded_queue_429_retry_after_and_slot_release(monkeypatch):
    """When active+queued jobs reach limit, new keys get 429 + Retry-After; slots release on success and error."""
    monkeypatch.setenv("JOB_MAX_ACTIVE", "2")
    gate = threading.Event()

    def blocking_exposures(ra, dec, det):
        gate.wait(timeout=5.0)
        if ra > 50:
            raise RuntimeError("boom")
        return []

    monkeypatch.setattr(irsa, "exposures", blocking_exposures)

    r1 = c.get("/api/view?ra=10.0&dec=10.0&fov=3&det=3&mode=coadd")
    r2 = c.get("/api/view?ra=60.0&dec=10.0&fov=3&det=3&mode=coadd")
    assert r1.status_code == 202 and r2.status_code == 202

    # 3rd distinct job while 2 are active -> must be rejected immediately with 429 or 503 and Retry-After
    r3 = c.get("/api/view?ra=30.0&dec=10.0&fov=3&det=3&mode=coadd")
    assert r3.status_code in (429, 503), f"Expected 429/503 when queue full, got {r3.status_code}"
    assert "Retry-After" in r3.headers
    assert "error" in r3.json()

    # Release the 2 running jobs (one succeeds, one raises RuntimeError)
    gate.set()
    for ra in (10.0, 60.0):
        for _ in range(30):
            if c.get(f"/api/view?ra={ra}&dec=10.0&fov=3&det=3&mode=coadd").status_code == 200:
                break
            time.sleep(0.02)

    # Now both slots must be free again
    r4 = c.get("/api/view?ra=30.0&dec=10.0&fov=3&det=3&mode=coadd")
    assert r4.status_code in (200, 202)


def test_partial_failure_not_cached_forever(monkeypatch):
    """If failed > 0 during cutout download, the incomplete result must not be locked in permanent cache forever."""
    monkeypatch.setenv("JOB_HOLD_SEC", "0.08")
    e0 = fake_exposure(CENTER, 100.0, seed=0)
    e1 = fake_exposure(CENTER, 101.0, seed=1)
    expos = {"u0": e0, "u1": e1}
    meta = [{"uri": "u0", "mjd": 100.0, "lam": [2.0, 2.1]}, {"uri": "u1", "mjd": 101.0, "lam": [2.0, 2.1]}]
    fail_second = {"on": True}

    def flaky_cutout(uri, ra, dec, size):
        if uri == "u1" and fail_second["on"]:
            raise ApiError(502, "upstream", "temporary cutout error")
        return uri

    monkeypatch.setattr(irsa, "exposures", lambda ra, dec, det: meta)
    monkeypatch.setattr(irsa, "cutout", flaky_cutout)
    monkeypatch.setattr(jobs, "load_exposure", lambda path, mjd_hint=None, obs_id_hint="": expos[path])

    url = f"/api/view?ra={CENTER.ra.deg}&dec={CENTER.dec.deg}&fov=3&det=3&mode=coadd"
    for _ in range(50):
        r = c.get(url)
        if r.status_code == 200:
            break
        time.sleep(0.01)
    d1 = r.json()
    assert d1["frames"][0]["n_exp"] == 1 and d1["warnings"]

    # Network recovers: next request after brief hold must re-run and get all 2 exposures
    fail_second["on"] = False
    time.sleep(0.12)
    for _ in range(50):
        r = c.get(url)
        if r.status_code == 200:
            break
        time.sleep(0.01)
    d2 = r.json()
    assert d2["frames"][0]["n_exp"] == 2 and not d2["warnings"]


def test_empty_result_cached_with_ttl_and_exp_list_ttl(monkeypatch, isolated_cache):
    """Empty coverage is cached for CACHE_EMPTY_TTL seconds (not dropped on first read), then refreshed."""
    monkeypatch.setenv("CACHE_EMPTY_TTL", "0.15")
    calls = {"n": 0}

    def exp_fn(ra, dec, det):
        calls["n"] += 1
        return []

    monkeypatch.setattr(irsa, "exposures", exp_fn)
    url = "/api/view?ra=88.0&dec=11.0&fov=3&det=3&mode=coadd"
    for _ in range(30):
        r = c.get(url)
        if r.status_code == 200:
            break
        time.sleep(0.02)
    assert calls["n"] == 1

    # Immediate second request must hit the TTL cache without calling irsa.exposures
    r2 = c.get(url)
    assert r2.status_code == 200 and calls["n"] == 1

    # After TTL expires, a new request queries irsa.exposures again
    time.sleep(0.2)
    for _ in range(30):
        r3 = c.get(url)
        if r3.status_code == 200:
            break
        time.sleep(0.02)
    assert calls["n"] == 2

    # Also test IRSA_EXP_TTL for _tap_exposures
    monkeypatch.setenv("IRSA_EXP_TTL", "0.1")
    tap_calls = {"n": 0}

    class TapResp:
        status_code = 200
        text = "uri,mjd,lo,hi\nlevel2_2025W10_1A_0001_1D3_spx.fits,60900.0,1e-6,2e-6\n"

    def fake_send(method, url, **kw):
        tap_calls["n"] += 1
        return TapResp()

    monkeypatch.setattr(irsa, "_send", fake_send)
    res1 = irsa._tap_exposures(77.0, 22.0, 3)
    res2 = irsa._tap_exposures(77.0, 22.0, 3)
    assert len(res1) == 1 and len(res2) == 1 and tap_calls["n"] == 1
    time.sleep(0.15)
    res3 = irsa._tap_exposures(77.0, 22.0, 3)
    assert len(res3) == 1 and tap_calls["n"] == 2


def test_parallel_spectrum_and_lightcurve_deduplication_and_upstream_error(monkeypatch, isolated_cache):
    """Parallel /api/spectrum and /api/lightcurve requests deduplicate to 1 run; upstream failure without stale returns error."""
    monkeypatch.setenv("JOB_HOLD_SEC", "0.5")
    calls = {"n": 0}
    barrier = threading.Barrier(4)

    def failing_exposures(ra, dec, det):
        if det == 1:
            calls["n"] += 1
        time.sleep(0.05)
        raise ApiError(502, "upstream", "IRSA down")

    monkeypatch.setattr(irsa, "exposures", failing_exposures)

    def poll_ep(path):
        barrier.wait()
        for _ in range(40):
            r = c.get(path)
            if r.status_code == 200:
                return r.json()
            time.sleep(0.02)
        raise AssertionError("timeout")

    paths = [
        "/api/spectrum?ra=19.0&dec=29.0&per_epoch=1",
        "/api/spectrum?ra=19.0&dec=29.0&per_epoch=1",
        "/api/lightcurve?ra=19.0&dec=29.0&per_epoch=1&det=2",
        "/api/lightcurve?ra=19.0&dec=29.0&per_epoch=1&det=2",
    ]
    with ThreadPoolExecutor(4) as ex:
        outs = [f.result() for f in [ex.submit(poll_ep, p) for p in paths]]

    assert all(o["status"] == "error" and "IRSA down" in o["message"] for o in outs)
    assert calls["n"] == 1, f"Expected 1 spectrum job run across spectrum+lightcurve, got {calls['n']}"


def test_upstream_env_backoff_and_stale_cache_fallback(monkeypatch, isolated_cache):
    """When archive fails with 502 upstream and stale cache exists, return stale cache with stale=True and warning."""
    monkeypatch.setenv("IRSA_BASE_URL", "http://127.0.0.1:9")
    monkeypatch.setenv("IRSA_RETRIES", "2")
    monkeypatch.setenv("IRSA_BACKOFF", "0.01")
    monkeypatch.setenv("SPHEREX_TABLES", str(isolated_cache / "no_tables"))

    # Seed an expired/old-version cached result for ra=45.0, dec=45.0, fov=5.0, det=3, mode=coadd, max_per_epoch=6
    k = cache.key(45.0, 45.0, 5.0, 3, "coadd", 6)
    d = isolated_cache / k
    d.mkdir(parents=True, exist_ok=True)
    (d / "result.json").write_text(json.dumps({
        "status": "ready", "query": {"ra": 45.0, "dec": 45.0}, "frames": [{"index": 0}], "warnings": [], "v": 1
    }))

    for _ in range(40):
        r = c.get("/api/view?ra=45.0&dec=45.0&fov=5&det=3&mode=coadd")
        if r.status_code == 200:
            break
        time.sleep(0.02)
    body = r.json()
    assert body["status"] == "ready"
    assert body.get("stale") is True
    assert any("устаревш" in w.lower() or "кэш" in w.lower() for w in body.get("warnings", []))


@pytest.mark.parametrize("q", [
    "ra=nan&dec=10",
    "ra=10&dec=inf",
    "ra=10&dec=-inf",
    "ra=10&dec=10&fov=nan",
    "ra=10&dec=10&fov=inf",
])
def test_nan_and_inf_rejected_with_400(q):
    r = c.get("/api/view?" + q)
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"


def test_lightcurve_nan_window_rejected():
    r = c.get("/api/lightcurve?ra=10&dec=10&wave_min=nan&wave_max=2.0")
    assert r.status_code == 400 and r.json()["error"]["code"] == "bad_param"
