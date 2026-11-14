"""Target presets (verbatim from the API spec), SIMBAD lookup via CDS Sesame, proper-motion propagation."""
from contextlib import contextmanager
import json
import logging
import os
import re
import threading
import time

import httpx
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from . import cache
from .errors import ApiError
from .params import Params

log = logging.getLogger("resolve")

# Section 3.2 of the spec, unchanged. No coordinates are stored here: they come from SIMBAD at request time.
PRESETS = [
    {"id": "barnard", "name": "Barnard's Star", "det": 3, "fov": 3, "note": "самое быстрое собственное движение, ~10.4″/год"},
    {"id": "epsindib", "name": "eps Ind B", "det": 3, "fov": 5, "note": "пара коричневых карликов, ищите метан"},
    {"id": "wise0855", "name": "WISE 0855-0714", "det": 6, "fov": 3, "note": "виден только на 4–5 мкм"},
    {"id": "rwaur", "name": "RW Aur", "det": 2, "fov": 3, "note": "переменная звезда T Tauri"},
    {"id": "rdra", "name": "R Dra", "det": 2, "fov": 3, "note": "мирида, период ~245 дней"},
]
# Sesame mirrors, tried in order (the primary host redirects to a flaky mirror).
SESAME = [
    "https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/SNV?",
    "https://vizier.idia.ac.za/viz-bin/nph-sesame/-oI/SNV?",
    "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?",
]
_mem = {}
_lock = threading.Lock()
_inflight: dict[str, threading.Lock] = {}
_FILE = cache.CACHE / "resolve.json"


@contextmanager
def _key_lock(k: str):
    with _lock:
        lk = _inflight.setdefault(k, threading.Lock())
    with lk:
        yield


def _sesame_urls() -> list[str]:
    env = os.environ.get("SESAME_URLS", "").strip()
    return [u.strip() for u in env.split(",") if u.strip()] if env else SESAME


def _fetch_sesame(name: str) -> list[str]:
    """Raw Sesame answers, one per mirror that responded. Raises 502 if none did."""
    out = []
    retries = int(os.environ.get("SESAME_RETRIES", "2"))
    backoff = float(os.environ.get("SESAME_BACKOFF", "0.5"))
    timeout = float(os.environ.get("SESAME_TIMEOUT", "30.0"))
    for base in _sesame_urls():
        for i in range(max(1, retries)):
            t = time.time()
            try:
                r = httpx.get(base + name, timeout=timeout, follow_redirects=True)
                log.info("GET %s -> %s %.1fs", base + name, r.status_code, time.time() - t)
                if r.status_code == 200:
                    out.append(r.text)
                    break
            except httpx.HTTPError as e:
                log.warning("sesame %s failed: %s", base, e)
            if i + 1 < retries:
                time.sleep(backoff * (2 ** i))
        if out and ("%J" in out[-1] or "Nothing found" in out[-1]):  # a definite answer: do not ask the other mirrors
            break
    if not out:
        raise ApiError(502, "upstream", "SIMBAD/Sesame недоступны")
    return out


def _simbad_section(text: str) -> list[str]:
    """Only the lines Sesame took from SIMBAD: NED/VizieR blocks may lack proper motion, never mix them in."""
    keep, on = [], False
    for line in text.splitlines():
        if line.startswith("#="):
            on = "=Sc=Simbad" in line
        elif on:
            keep.append(line)
    return keep


def _parse(name: str, texts: list[str]) -> dict | None:
    for text in texts:
        lines = _simbad_section(text)
        j = next((re.match(r"%J\s+([\d.+-]+)\s+([\d.+-]+)", l) for l in lines if l.startswith("%J ")), None)
        if not j:
            continue
        pm = next((re.match(r"%P\s+([\d.+-]+)\s+([\d.+-]+)", l) for l in lines if l.startswith("%P ")), None)
        x = next((re.match(r"%X\s+([\d.+-]+)", l) for l in lines if l.startswith("%X ")), None)
        return {"name": name, "ra": float(j[1]), "dec": float(j[2]),
                "pmra_masyr": float(pm[1]) if pm else None, "pmdec_masyr": float(pm[2]) if pm else None,
                "parallax_mas": float(x[1]) if x else None, "ref_epoch": 2000.0,
                "source": "SIMBAD (via CDS Sesame)"}
    return None


_team = {}  # OBJECT_PARAMETERS.csv: mtime and parsed targets


def _preset_record(query: str) -> dict | None:
    """The five challenge targets come from the team's OBJECT_PARAMETERS.csv (Gaia DR3 / Kirkpatrick+2021).
    That table already settles what SIMBAD names get wrong: WISE 0855 (the name is not a SIMBAD identifier),
    RW Aur A/B, eps Ind Ba. Matching by id or name is done by src/processing/targets.find_preset."""
    from processing.targets import DEFAULT_CATALOG, find_preset, load_targets

    if not DEFAULT_CATALOG.exists():
        return None
    m = DEFAULT_CATALOG.stat().st_mtime
    if _team.get("m") != m:
        _team.update(m=m, targets=load_targets())
    t = find_preset(query, _team["targets"])
    if t is None:
        return None
    return {"name": t.name, "id": t.id, "ra": t.ra, "dec": t.dec, "pmra_masyr": t.pmra, "pmdec_masyr": t.pmdec,
            "parallax_mas": t.parallax or None, "rv_kms": t.rv or None, "ref_epoch": float(t.ref_epoch.jyear),
            "source": t.source}


def resolve(name: str) -> dict:
    """Preset id or name -> team table; any other name -> SIMBAD (cached forever: memory + resolve.json).
    Result: position at ref_epoch, proper motion, parallax."""
    rec = _preset_record(name.strip())
    if rec:
        return rec
    k = name.strip().lower()
    with _key_lock(k):
        with _lock:
            if not _mem and _FILE.exists():
                try:
                    _mem.update(json.loads(_FILE.read_text(encoding="utf-8")))
                except Exception:
                    pass
            if k in _mem:
                return _mem[k]
        rec = _parse(name.strip(), _fetch_sesame(name.strip()))
        if rec is None:
            raise ApiError(404, "not_found", f"объект {name!r} не найден в SIMBAD")
        with _lock:
            _mem[k] = rec
            cache.atomic_write_text(_FILE, json.dumps(_mem, ensure_ascii=False))
        return rec


def propagate(o: dict, mjd: float):
    """Position at `mjd` from J2000 position + proper motion (+ parallax distance). No motion known: unchanged."""
    if o.get("pmra_masyr") is None or o.get("pmdec_masyr") is None:
        return o["ra"], o["dec"]
    kw = dict(ra=o["ra"] * u.deg, dec=o["dec"] * u.deg,
              pm_ra_cosdec=o["pmra_masyr"] * u.mas / u.yr, pm_dec=o["pmdec_masyr"] * u.mas / u.yr,
              obstime=Time(o["ref_epoch"], format="jyear"))
    plx = o.get("parallax_mas")
    if plx and plx > 0:
        kw["distance"] = (1000.0 / plx) * u.pc
        kw["radial_velocity"] = (o.get("rv_kms") or 0.0) * u.km / u.s
    else:
        kw["radial_velocity"] = 0 * u.km / u.s
    c = SkyCoord(**kw).apply_space_motion(new_obstime=Time(mjd, format="mjd"))
    return float(c.ra.deg), float(c.dec.deg)


def center_mjd() -> float:
    """Date the query centre is propagated to: today. CENTER_MJD pins it (keeps cache keys stable for a demo)."""
    return float(os.environ["CENTER_MJD"]) if os.environ.get("CENTER_MJD") else float(Time.now().mjd)


def for_params(p: Params) -> dict:
    """-> {center:(ra,dec) used for the IRSA query and the grid, motion: SIMBAD record or None, label}."""
    if p.ra is not None:  # explicit coordinates: no proper motion is known, so no predicted position either
        return {"center": (p.ra, p.dec), "motion": None, "label": p.target}
    o = resolve(p.target)
    return {"center": propagate(o, center_mjd()), "motion": o, "label": o["name"]}
