"""Asynchronous API jobs for on-the-fly SPHEREx spectra and light curves."""
from __future__ import annotations

import hashlib
import json
import logging
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord

from . import cache, resolve
from .errors import ApiError
from .params import parse_coords

log = logging.getLogger("science_jobs")
# The project keeps reusable science code in src/; make it importable when
# uvicorn is launched from the repository root (as in start.sh).
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from processing.pipeline import build_spectrum, lightcurve_from_spectrum  # noqa: E402

_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="spherex-science")
_LOCK = threading.Lock()
_JOBS: dict[str, dict] = {}


def _target(target: str | None, ra: str | None, dec: str | None):
    if (ra is None) != (dec is None):
        raise ApiError(400, "bad_param", "ra и dec нужно передавать вместе")
    if ra is not None:
        r, d = parse_coords(ra, dec)
        coord = SkyCoord(r * u.deg, d * u.deg)
        return coord, (lambda _mjd: coord), f"{r:.7f},{d:.7f}", None
    if not target:
        raise ApiError(400, "bad_param", "нужен target либо пара ra+dec")
    # Accept both the preset IDs used in shareable URLs (e.g. "barnard")
    # and arbitrary object names entered by the user.
    preset = next((item for item in resolve.PRESETS
                   if item["id"] == target.strip().lower()), None)
    lookup_name = preset["name"] if preset else target.strip()
    obj = resolve.resolve(lookup_name)
    now = resolve.center_mjd()
    r, d = resolve.propagate(obj, now)
    center = SkyCoord(r * u.deg, d * u.deg)

    def coord_at(mjd: float):
        rr, dd = resolve.propagate(obj, mjd)
        return SkyCoord(rr * u.deg, dd * u.deg)

    key = f"{lookup_name.strip().lower()}|{obj.get('ra')}|{obj.get('dec')}|{obj.get('pmra_masyr')}|{obj.get('pmdec_masyr')}"
    return center, coord_at, key, obj


def _key(target_key: str) -> str:
    return hashlib.sha1(target_key.encode("utf-8")).hexdigest()[:16]


def _path(key: str) -> Path:
    return cache.CACHE / f"spectrum_{key}.json"


def _run(key: str, center: SkyCoord, coord_at, target_info):
    try:
        result = build_spectrum(center, coord_at)
        result["target"] = target_info
        tmp = _path(key).with_suffix(".json.tmp")
        tmp.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False))
        tmp.replace(_path(key))
        with _LOCK:
            _JOBS[key] = {"status": "ready", "result": result}
    except Exception as exc:  # keep error visible to the polling client
        log.exception("spectrum job %s failed", key)
        with _LOCK:
            _JOBS[key] = {"status": "error", "message": str(exc)}


def spectrum(target: str | None, ra: str | None, dec: str | None):
    center, coord_at, target_key, target_info = _target(target, ra, dec)
    key = _key(target_key)
    path = _path(key)
    if path.exists():
        return 200, json.loads(path.read_text())
    with _LOCK:
        job = _JOBS.get(key)
        if job is None:
            _JOBS[key] = {"status": "running"}
            _POOL.submit(_run, key, center, coord_at, target_info)
            return 202, {"status": "running", "job_id": key, "message": "Считаем спектр SPHEREx"}
        if job["status"] == "running":
            return 202, {"status": "running", "job_id": key, "message": "Считаем спектр SPHEREx"}
        if job["status"] == "error":
            message = job["message"]
            _JOBS.pop(key, None)
            raise ApiError(502, "processing", message)
        return 200, job["result"]


def lightcurve(target: str | None, ra: str | None, dec: str | None,
               wave_min: float, wave_max: float):
    if not (0 < wave_min < wave_max <= 5.1):
        raise ApiError(400, "bad_param", "нужен диапазон 0 < wave_min < wave_max <= 5.1 мкм")
    status, spec = spectrum(target, ra, dec)
    if status == 202:
        return status, {**spec, "message": "Сначала строится спектр; повторите запрос для кривой блеска"}
    result = lightcurve_from_spectrum(spec, wave_min, wave_max)
    result["target"] = spec.get("target")
    return 200, result
 
