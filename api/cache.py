"""Result cache: CACHE_DIR/<key>/result.json + PNGs. FITS cutouts and lookups live in CACHE_DIR too."""
import hashlib
import json
import os
from pathlib import Path

CACHE = Path(os.environ.get("CACHE_DIR", Path(__file__).resolve().parent.parent / "cache"))
CACHE.mkdir(parents=True, exist_ok=True)


def key(ra, dec, fov, det, mode, max_per_epoch) -> str:
    s = f"{round(ra, 5)}|{round(dec, 5)}|{fov}|{det}|{mode}|{max_per_epoch}"
    return hashlib.sha1(s.encode()).hexdigest()[:16]


VERSION = 4  # bump when the result format changes: older cached results are then recomputed (FITS stay cached)


def key_spectrum(ra, dec) -> str:
    return hashlib.sha1(f"spectrum|{round(ra, 5)}|{round(dec, 5)}".encode()).hexdigest()[:16]


def read(k: str):
    p = CACHE / k / "result.json"
    if not p.exists():
        return None
    res = json.loads(p.read_text(encoding="utf-8"))
    return res if res.pop("v", None) == VERSION else None  # "v" is internal, never sent to clients


def write(k: str, obj: dict):
    """result.json is written last and atomically: its presence means the PNGs are complete."""
    d = CACHE / k
    d.mkdir(exist_ok=True)
    tmp = d / "result.json.tmp"
    tmp.write_text(json.dumps({**obj, "v": VERSION}, ensure_ascii=False), encoding="utf-8")
    tmp.replace(d / "result.json")
