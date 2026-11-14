"""Result cache: CACHE_DIR/<key>/result.json + PNGs. FITS cutouts and lookups live in CACHE_DIR too."""
import hashlib
import json
import os
import time
import uuid
from pathlib import Path

CACHE = Path(os.environ.get("CACHE_DIR", Path(__file__).resolve().parent.parent / "cache"))
CACHE.mkdir(parents=True, exist_ok=True)


def key(ra, dec, fov, det, mode, max_per_epoch) -> str:
    s = f"{round(ra, 5)}|{round(dec, 5)}|{fov}|{det}|{mode}|{max_per_epoch}"
    return hashlib.sha1(s.encode()).hexdigest()[:16]


VERSION = 4  # bump when the result format changes: older cached results are then recomputed (FITS stay cached)


def atomic_write_bytes(path: Path, data: bytes):
    """Write bytes atomically via a unique temporary file in the same directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f".tmp_{uuid.uuid4().hex[:8]}_{path.name}"
    tmp.write_bytes(data)
    tmp.replace(path)


def atomic_write_text(path: Path, text: str):
    """Write text atomically via a unique temporary file in the same directory."""
    atomic_write_bytes(path, text.encode("utf-8"))


def read(k: str, allow_stale: bool = False):
    p = CACHE / k / "result.json"
    if not p.exists():
        return None
    try:
        res = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None
    v = res.pop("v", None)
    exp = res.pop("expires_at", None)
    if allow_stale:
        return res
    if v != VERSION:
        return None
    if exp is not None and time.time() >= exp:
        return None
    return res


def write(k: str, obj: dict, ttl: float | None = None):
    """result.json is written last and atomically: its presence means the PNGs are complete."""
    d = CACHE / k
    d.mkdir(parents=True, exist_ok=True)
    payload = {**obj, "v": VERSION}
    if ttl is not None:
        payload["expires_at"] = time.time() + float(ttl)
    atomic_write_text(d / "result.json", json.dumps(payload, ensure_ascii=False))
