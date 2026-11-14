"""Fill the API cache for the demo by calling the running API: all presets, the shown modes, spectrum, light curve.

The server MUST run with the same CACHE_DIR and CENTER_MJD that will be used at the demo (see start_demo.sh):
the cache key depends on the query centre, and the centre is propagated to CENTER_MJD (default: today).

    bash start_demo.sh &                                   # terminal 1
    .venv/bin/python scripts/prepare_demo_cache.py         # terminal 2
Options: --base URL, --modes coadd,diff,trail,single, --only barnard,rwaur, --retries N, --report file.json
"""
import argparse
import json
import sys
import threading
import time

import httpx


def poll(base, path, params, limit=1800):
    """GET until 200 (202 = job running). -> (status, body, seconds)."""
    t0 = time.time()
    while time.time() - t0 < limit:
        r = httpx.get(base + path, params=params, timeout=300)
        body = r.json()
        if r.status_code == 202:
            time.sleep(5)
            continue
        return r.status_code, body, round(time.time() - t0)
    return 0, {"status": "error", "message": "timeout"}, limit


def run_one(base, path, params, retries):
    """IRSA answers 503 on about half of the requests: a failed job is simply asked again."""
    for attempt in range(1, retries + 1):
        status, body, secs = poll(base, path, params)
        ok = status == 200 and body.get("status") == "ready" and not body.get("error")
        if ok:
            return {"ok": True, "secs": secs, "attempt": attempt, "body": body}
        print(f"  retry {attempt}/{retries}: {path} {params} -> {status} {body.get('message') or body.get('error')}",
              flush=True)
    return {"ok": False, "secs": secs, "attempt": retries, "body": body}


def summarize(path, body):
    if path == "/api/view":
        f = body.get("frames", [])
        return f"{len(f)} frames, trail={'yes' if body.get('trail_url') else 'no'}, warnings={len(body.get('warnings', []))}"
    if path == "/api/spectrum":
        pts = body.get("points", [])
        usable = sum(1 for p in pts if p.get("flux_ujy") is not None)
        return f"{len(pts)} points, {usable} with flux, {body.get('n_rejected')} flagged"
    return f"{len(body.get('points', []))} points in {body.get('wave_window_um')}"


def target_worker(base, t, modes, retries, results, lock):
    jobs = [("/api/view", {"target": t["id"], "fov": t["fov"], "det": t["det"], "mode": m}) for m in modes]
    jobs += [("/api/spectrum", {"target": t["id"]}), ("/api/lightcurve", {"target": t["id"], "det": t["det"]})]
    for path, params in jobs:
        res = run_one(base, path, params, retries)
        line = f"{t['id']:9s} {path:15s} {params.get('mode', ''):6s} "
        line += (f"OK    {res['secs']:4d}s  " + summarize(path, res["body"])) if res["ok"] else "FAILED"
        print(line, flush=True)
        with lock:
            results.append({"target": t["id"], "path": path, "params": params, "ok": res["ok"], "secs": res["secs"],
                            "summary": summarize(path, res["body"]) if res["ok"] else str(res["body"])})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://127.0.0.1:8000")
    ap.add_argument("--modes", default="coadd,diff")
    ap.add_argument("--only", default="")
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--report", default="demo_cache_report.json")
    a = ap.parse_args()
    targets = httpx.get(a.base + "/api/targets", timeout=60).json()
    if a.only:
        targets = [t for t in targets if t["id"] in a.only.split(",")]
    results, lock = [], threading.Lock()
    threads = [threading.Thread(target=target_worker, args=(a.base, t, a.modes.split(","), a.retries, results, lock))
               for t in targets]
    [t.start() for t in threads]
    [t.join() for t in threads]
    json.dump(results, open(a.report, "w"), ensure_ascii=False, indent=2)
    bad = [r for r in results if not r["ok"]]
    print(f"\n{len(results) - len(bad)} of {len(results)} requests ready; report: {a.report}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
