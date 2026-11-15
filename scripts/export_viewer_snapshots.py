"""Build every viewer mode from the exact saved real FITS, with network disabled."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "scripts")]
os.environ["PREWARM"] = "0"

from astropy.utils import iers
from api import cache, irsa, jobs, params, resolve
from export_challenge_gifs import SNAPSHOTS, dump, no_network, saved_input

OUT = ROOT / "frontend/challenge-demo/snapshots"


def export_query(query, rows, lookup, manifest):
    p = params.parse(**query)
    resolved = resolve.for_params(p, epoch=61358)
    key = jobs._view_key(p, resolved, source="archive")
    metadata = irsa.ExposureRows(rows, source="saved_sia2_replay", warnings=[
        "Сохранённая выборка IRSA, повторно обработанная локально; новый поиск архива не выполнялся"])
    with patch.object(irsa, "exposures", lambda *a: metadata), \
         patch.object(irsa, "cutout", lambda uri, *a: lookup[uri]):
        result = jobs._pipeline(jobs.Job(key, "snapshot"), p, resolved)
    assert result["frames"] and result["exposure_counts"]["failed"] == 0, query
    dest = OUT / key
    dest.mkdir(parents=True, exist_ok=True)
    urls = [(frame, field) for frame in result["frames"] for field in ("image_url", "flags_url")]
    if result.get("trail_url"):
        urls.append((result, "trail_url"))
    for record, field in urls:
        source = cache.CACHE / record[field].removeprefix("/cache/")
        shutil.copyfile(source, dest / source.name)
        record[field] = f"/challenge-demo/snapshots/{key}/{source.name}"
    result.update(snapshot_date="2026-11-15", snapshot_version=cache.VERSION)
    dump(dest / "result.json", result)
    entry = {"key": key, "request": query, "frames": len(result["frames"]),
             "exposures": result["exposure_counts"]["used"], "inputs_manifest": manifest,
             "result_sha256": hashlib.sha256((dest / "result.json").read_bytes()).hexdigest()}
    print(json.dumps({k: entry[k] for k in ("request", "frames", "exposures")}), flush=True)
    return entry


def run():
    iers.conf.auto_download = False
    entries = []
    with tempfile.TemporaryDirectory(prefix="spherex-snapshot-export-") as temporary, \
         patch.object(cache, "CACHE", Path(temporary)), \
         patch.object(irsa, "_send", no_network), patch("httpx.get", no_network), \
         patch("httpx.request", no_network), patch("requests.sessions.Session.request", no_network), \
         patch("socket.create_connection", no_network):
        for target in SNAPSHOTS:
            old, rows, lookup, inputs, source = saved_input(target)
            original = json.loads((ROOT / f"frontend/challenge-demo/{target}/inputs.json").read_text())
            assert {(i["product_url"], i["sha256"]) for i in inputs} == {
                (i["product_url"], i["sha256"]) for i in original["inputs"]}
            for maximum in (6, 1):
                for mode in params.MODES:
                    entries.append(export_query({"target": target, "fov": old["query"]["fov"],
                        "det": old["query"]["det"], "mode": mode, "max_per_epoch": maximum},
                        rows, lookup, f"{target}/inputs.json"))
        # Preserve the exact coordinate example as well as the five tracking presets.
        manifest = json.loads((ROOT / "frontend/challenge-demo/wise0855-coordinates/inputs.json").read_text())
        inputs = manifest["inputs"]
        lookup = {i["product_url"]: ROOT / i["path"] for i in inputs}
        assert all(hashlib.sha256(lookup[i["product_url"]].read_bytes()).hexdigest() == i["sha256"] for i in inputs)
        rows = [{"uri": i["product_url"], "mjd": i["mjd_utc"]} for i in inputs]
        for maximum in (6, 1):
            for mode in params.MODES:
                entries.append(export_query({"ra": "08h55m10.831s", "dec": "-7d14m42.53s", "fov": 6,
                    "det": 6, "mode": mode, "max_per_epoch": maximum}, rows, lookup,
                    "wise0855-coordinates/inputs.json"))
    dump(OUT / "manifest.json", {"snapshot_date": "2026-11-15", "cache_version": cache.VERSION,
        "center_mjd": 61358, "network_requests": 0, "entries": entries,
        "processing_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                              for directory in ("api", "src/processing") for p in (ROOT / directory).glob("*.py")}})


if __name__ == "__main__":
    run()
