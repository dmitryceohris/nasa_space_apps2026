"""Result cache: CACHE_DIR/<key>/result.json + PNGs. FITS cutouts and lookups live in CACHE_DIR too."""
import hashlib
import json
import os
import time
from pathlib import Path

CACHE = Path(os.environ.get("CACHE_DIR", Path(__file__).resolve().parent.parent / "cache"))
CACHE.mkdir(parents=True, exist_ok=True)


def key(ra, dec, fov, det, mode, max_per_epoch) -> str:
    s = f"{round(ra, 5)}|{round(dec, 5)}|{fov}|{det}|{mode}|{max_per_epoch}"
    return hashlib.sha1(s.encode()).hexdigest()[:16]


VERSION = 12  # Scientific validation and corrected masks; old measurements must be recomputed.


def key_spectrum(ra, dec) -> str:
    return hashlib.sha1(f"spectrum|{round(ra, 5)}|{round(dec, 5)}".encode()).hexdigest()[:16]


def read(k: str):
    p = CACHE / k / "result.json"
    try:
        res = json.loads(p.read_text(encoding="utf-8"))
        modified = p.stat().st_mtime
    except (OSError, ValueError):
        return None
    if not isinstance(res, dict) or res.pop("v", None) != VERSION:
        return None
    return res if time.time() - modified < ttl(res) else None


def ttl(res):
    empty = ("frames" in res and not res["frames"]) or ("points" in res and not res["points"])
    # Any state that may come from failed downloads is retried soon (spectrum: unusable_data can keep points).
    retry = res.get("data_status") in ("partial", "download_failed", "unusable_data")
    setting = "EMPTY_CACHE_TTL" if empty else "PARTIAL_CACHE_TTL" if retry else "RESULT_CACHE_TTL"
    return float(os.environ.get(setting, "60" if setting != "RESULT_CACHE_TTL" else "3600"))


def write(k: str, obj: dict):
    """result.json is written last and atomically: its presence means the PNGs are complete."""
    d = CACHE / k
    d.mkdir(exist_ok=True)
    tmp = d / "result.json.tmp"
    tmp.write_text(json.dumps({**obj, "v": VERSION}, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    tmp.replace(d / "result.json")
