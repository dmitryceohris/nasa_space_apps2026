"""Replay saved real cutouts through /api/view without archive traffic or shared-cache writes.

Only archive selection/download and Sesame transport are replaced by saved inputs.
FITS loading, masks, WCS, background, coadd/diff/trail, jobs, JSON and PNG serving are real.
This checks a bounded historical sample, not archive completeness or source detection.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
os.environ["PREWARM"] = "0"

import httpx
from astropy.utils import iers
from fastapi.testclient import TestClient

from Exposure import load_exposure
from api import cache, irsa, resolve, main
from api.limits import RequestLimiter
from api.main import app
from processing.frames import lone_spike_mask


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def no_network(*args, **kwargs):
    raise AssertionError("The acceptance replay must not access the archive/network")


def run(fits_cache, download_limit=0):
    logging.getLogger("httpx").setLevel(logging.WARNING)
    evidence = json.loads((ROOT / "data/validation/science-acceptance/evidence.json").read_text())
    rows = evidence["frames"] + evidence["wise_exact_coordinate_frames"]
    needed = {r["cutout_sha256"] for r in rows}
    paths = {}
    for path in fits_cache.glob("cut_*.fits"):
        digest = sha(path)
        if digest in needed:
            paths[digest] = path
    missing = needed - paths.keys()
    assert not missing or len(missing) <= download_limit, (
        f"Missing {len(missing)} exact saved cutouts. Supply the saved FITS directory or explicitly "
        f"allow at most {len(missing)} single-attempt downloads with --download-limit {len(missing)}")
    # Optional reproduction on another machine: exact manifest URLs only, no archive search/retries.
    if missing:
        fits_cache.mkdir(parents=True, exist_ok=True)
        for digest in sorted(missing):
            row = next(r for r in rows if r["cutout_sha256"] == digest)
            response = httpx.get(row["cutout_url"], timeout=45, follow_redirects=True)
            response.raise_for_status()
            assert hashlib.sha256(response.content).hexdigest() == digest, row["obs_id"]
            path = fits_cache / f"cut_{digest}.fits"
            path.write_bytes(response.content)
            paths[digest] = path
    assert needed <= paths.keys(), f"Missing {len(needed - paths.keys())} exact saved cutouts"
    for source in json.loads((ROOT / "data/validation/science-acceptance/sources.json").read_text()):
        assert sha(ROOT / "data/validation/science-acceptance" / source["file"]) == source["sha256"]
    products = list(csv.DictReader((ROOT / "data/validation/science-acceptance/irsa-products.csv").open()))
    products = {r["uri"]: r for r in products}
    input_checks = []
    for row in rows:
        assert sha(ROOT / row["response"]) == row["response_sha256"]
        exp = load_exposure(str(paths[row["cutout_sha256"]]))
        assert exp.obs_id == row["obs_id"] and exp.detector == row["detector"]
        assert exp.time_source == "MJD-AVG" and exp.time_assumptions == []
        assert exp.mjd == row["fits_header"]["MJD-AVG"]
        assert exp.header["DATE-AVG"] == row["date_utc"]
        uri = row["product_url"].removeprefix(irsa.IRSA + "/")
        assert irsa.obs_id(uri) == exp.obs_id
        assert products[uri]["energy_bandpassname"] == f"SPHEREx-D{exp.detector}"
        assert abs(float(products[uri]["time_bounds_lower"]) - exp.header["MJD-BEG"]) < 1e-6
        assert f"/spherex/{row['release']}/" in row["product_url"]
        x, y = (int(round(v)) for v in row["native_xy"])
        flags = exp.flags[y-2:y+3, x-2:x+3]
        assert flags.tolist() == row["flags_5x5_values"]
        spikes = lone_spike_mask(exp)
        assert not (spikes & exp.flag_mask(("SOURCE",))).any()
        for name, count in row["flag_counts_5x5"].items():
            assert int(((flags & (1 << exp.flag_bits[name])) != 0).sum()) == count
        input_checks.append({"response": row["response"], "obs_id": exp.obs_id,
                             "det": exp.detector, "release": row["release"],
                             "sha256": row["cutout_sha256"], "mjd_utc": exp.mjd,
                             "native_spikes_rejected": int(spikes.sum()),
                             "flags_5x5": row["flag_counts_5x5"]})

    results = []
    iers.conf.auto_download = False
    with tempfile.TemporaryDirectory(prefix="spherex-acceptance-") as temporary:
        output = Path(temporary)
        mount = next(r.app for r in app.routes if getattr(r, "name", "") == "cache")
        with patch.object(cache, "CACHE", output), patch.object(resolve, "_FILE", output / "resolve.json"), \
             patch.object(resolve, "_mem", {}), patch.object(mount, "all_directories", [str(output)]), \
             patch.object(main, "request_limiter", RequestLimiter(100000)), \
             patch.object(irsa, "_send", no_network), patch("httpx.request", no_network):
            with TestClient(app) as client:
                for response_path in dict.fromkeys(r["response"] for r in rows):
                    saved = json.loads((ROOT / response_path).read_text())
                    exact = "wise-objective" in response_path
                    old = saved["response"] if exact else saved["responses"]["view"]
                    saved_input = "WISE objective coordinates" if exact else saved["input"]
                    selected = [r for r in rows if r["response"] == response_path]
                    lookup = {r["product_url"]: paths[r["cutout_sha256"]] for r in selected}
                    metadata = [{"uri": r["product_url"], "mjd": float(products[
                        r["product_url"].removeprefix(irsa.IRSA + "/")]["time_bounds_lower"])} for r in selected]
                    query = {"fov": 6, "det": selected[0]["detector"], "max_per_epoch": 1}
                    if exact:
                        query.update(ra="08h55m10.831s", dec="-7d14m42.53s")
                    else:
                        query["target"] = saved_input
                    sesame_path = {"Proxima Cen": "proxima-cen", "61 Cyg A": "61-cyg-a", "Vega": "vega"}.get(saved_input)
                    texts = [(ROOT / f"data/validation/responses/{sesame_path}-sesame.txt").read_text()] if sesame_path else []
                    with patch.dict(os.environ, {"CENTER_MJD": str(old.get("center_mjd") or 61358),
                                                 "API_RATE_LIMIT_PER_MINUTE": "10000"}), \
                         patch.object(irsa, "exposures", lambda *args: metadata), \
                         patch.object(irsa, "cutout", lambda uri, *args: lookup[uri]), \
                         patch.object(resolve, "_fetch_sesame", lambda name: texts):
                        if not exact:
                            actual = client.get("/api/resolve", params={"target": saved_input}).json()
                            expected = saved["responses"]["resolve"]
                            for key in ("source_id", "ref_epoch", "ref_epoch_format", "ref_epoch_scale", "coordinate_model"):
                                assert actual[key] == expected[key], (saved_input, key)
                        for mode in ("single", "coadd", "diff", "trail"):
                            query["mode"] = mode
                            started = time.monotonic()
                            polls = 0
                            while True:
                                response = client.get("/api/view", params=query)
                                polls += 1
                                if response.status_code != 202:
                                    break
                                assert time.monotonic() - started < 60
                                time.sleep(.01)
                            assert response.status_code == 200, response.text
                            result = response.json()
                            assert result["status"] == "ready", result
                            assert len(result["frames"]) == len(selected)
                            assert result["exposure_counts"]["used"] == len(selected)
                            assert not result["processing"]["registration"]["independently_verified"]
                            if exact:
                                assert result["target_info"] is None
                            else:
                                assert result["motion"]["pm_measured_arcsec_yr"] is None
                            by_id = {r["obs_id"]: r for r in selected}
                            for frame in result["frames"]:
                                assert len(frame["obs_ids"]) == 1
                                ref = by_id[frame["obs_ids"][0]]
                                assert abs(frame["mjd"] - ref["fits_header"]["MJD-AVG"]) <= .00005
                                assert frame["date"] == ref["date_utc"][:10]
                                assert abs(frame["wave_um"] - ref["wave_um"]) <= .000051
                                assert frame["quality"]["verified"] is False
                                assert frame["time_sources"] == [{"obs_id": ref["obs_id"], "source": "MJD-AVG", "assumptions": []}]
                                for field in ("image_url", "flags_url"):
                                    image = client.get(frame[field])
                                    assert image.status_code == 200 and image.content.startswith(b"\x89PNG")
                                if mode == "diff":
                                    assert all(c["sign"] in (-1, 1) for c in frame["measured_candidates"])
                                if frame["measurement"]["target_association"] == "far_from_prediction":
                                    assert frame["measured_px"] is None
                            if mode in ("diff", "trail"):
                                assert result["processing"]["difference"]["units"] == "robust_residual_scale"
                                image = client.get(result["trail_url"])
                                assert image.status_code == 200 and image.content.startswith(b"\x89PNG")
                            results.append({"input": saved_input, "exact_wise_coordinates": exact,
                                            "mode": mode, "polls": polls, "response": result})
    return {"checked_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "scope": __doc__, "archive_requests": len(missing), "checked_inputs": input_checks,
            "preserved_inputs_sha256": {p: sha(ROOT / p) for p in (
                "OBJECT_PARAMETERS.csv", "output/predictions.csv", "data/physics_observations.csv")},
            "scenarios": results, "summary": {"real_inputs": len(input_checks),
                "http_scenarios": len(results), "rendered_frames": sum(len(r["response"]["frames"]) for r in results),
                "single_non_null_candidate_markers": sum(f["measured_px"] is not None for r in results
                    if r["mode"] == "single" for f in r["response"]["frames"]),
                "single_measurement_status": dict(Counter(f["measurement"]["status"] for r in results
                    if r["mode"] == "single" for f in r["response"]["frames"]))}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits-cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--download-limit", default=0, type=int, choices=range(48), metavar="0..47",
                        help="Explicit cap on exact-URL downloads of missing cutouts; default 0 (offline)")
    args = parser.parse_args()
    result = run(args.fits_cache, args.download_limit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps(result["summary"], ensure_ascii=False))
