"""Replay public-final view requests through GET /api/view with the SAME obs_ids, served from the local FITS
cache (IRSA TAP answered 502). Only the archive lookup/download is replaced; resolve, processing and the HTTP
path are real. Run from the repository root: uv run --locked python <this> <out.json>"""
import glob
import json
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path

os.environ["PREWARM"] = "0"
os.environ["CENTER_MJD"] = "61358"
sys.path[:0] = [".", "src"]
warnings.filterwarnings("ignore")

from fastapi.testclient import TestClient  # noqa: E402

from api import cache, irsa  # noqa: E402
from api.main import app  # noqa: E402
from Exposure import load_exposure  # noqa: E402

cache.CACHE = Path(tempfile.mkdtemp(prefix="replay_"))  # results only; FITS come from the project cache
by_obs = {}
for p in sorted(glob.glob("cache/cut_*.fits")):
    try:
        e = load_exposure(p)
    except Exception:  # noqa: BLE001
        continue
    key = (e.obs_id, e.detector)
    if key not in by_obs or max(e.shape) > by_obs[key][1]:
        by_obs[key] = (p, max(e.shape), e.mjd)

client = TestClient(app)
out = {}
for fname, name in (("0-barnard.json", "barnard"), ("1-epsindib.json", "epsindib"), ("2-wise0855.json", "wise0855")):
    stored = json.load(open(f"data/validation/public-final/responses/{fname}", encoding="utf-8"))["responses"]["view"]
    det = stored["query"]["det"]
    obs = [fr["obs_ids"][0] for fr in stored["frames"]]
    have = [o for o in obs if (o, det) in by_obs]
    meta = [{"uri": by_obs[(o, det)][0], "mjd": by_obs[(o, det)][2], "lam": [0, 0]} for o in have]
    irsa.exposures = lambda ra, dec, d, meta=meta: meta
    irsa.cutout = lambda uri, *a: uri
    q = {"target": name, "fov": stored["query"]["fov"], "det": det, "mode": "single", "max_per_epoch": 1}
    for _ in range(600):
        r = client.get("/api/view", params=q)
        if r.status_code != 202:
            break
        time.sleep(0.2)
    res = r.json()
    old = {fr["obs_ids"][0]: fr for fr in stored["frames"]}
    rows = []
    for fr in res.get("frames", []):
        o = old.get(fr["obs_ids"][0], {})
        m = fr["measurement"]
        rows.append({"obs_id": fr["obs_ids"][0], "date": fr["date"], "wave_um": fr["wave_um"],
                     "coverage": fr["quality"]["coverage_fraction"],
                     "old_measured_px": o.get("measured_px"), "new_measured_px": fr["measured_px"],
                     "predicted_px": fr["predicted_px"], "status": m["status"], "quality": m["quality"],
                     "association": m.get("target_association"), "offset_px": m.get("offset_from_prediction_px"),
                     "flags": m["flags"], "candidates": [(c["kind"], c["sign"], c["xy"], c["quality"])
                                                         for c in fr["measured_candidates"]]})
    out[name] = {"stored_obs": obs, "replayed_obs": have, "http": r.status_code,
                 "spectral_sampling": res.get("spectral_sampling"), "warnings": res.get("warnings"), "frames": rows}
    print(f"== {name}: replayed {len(have)}/{len(obs)} obs, HTTP {r.status_code}", flush=True)
    print("   spectral_sampling", res.get("spectral_sampling"), flush=True)
    for w in res.get("warnings", []):
        print("   warning:", w[:160], flush=True)
    for x in rows:
        print(f"   {x['date']} {x['obs_id']} λ={x['wave_um']} cov={x['coverage']:.2f} old={x['old_measured_px']} "
              f"new={x['new_measured_px']} pred={x['predicted_px']} {x['status']}/{x['quality']} "
              f"{x['association']} off={x['offset_px']} {x['flags']}", flush=True)
Path(sys.argv[1]).write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
