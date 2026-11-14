"""Catalogue-first target identity, Sesame fallback and the common astrometric adapter."""
import json
import logging
import math
import os
import re
import threading
import time
from urllib.parse import quote

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
    {"id": "wise0855", "name": "WISE 0855-0714", "det": 6, "fov": 3, "note": "D6: 4.42–5.00 мкм; обнаружение требует проверки"},
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
_FILE = cache.CACHE / "resolve.json"


def _fetch_sesame(name: str) -> list[str]:
    """Raw Sesame answers, one per mirror that responded. Raises 502 if none did."""
    out = []
    for base in SESAME:
        t = time.time()
        try:
            r = httpx.get(base + quote(name, safe=""), timeout=30, follow_redirects=True)
            log.info("GET %s -> %s %.1fs", base + quote(name, safe=""), r.status_code, time.time() - t)
            if r.status_code == 200:
                if "%J" in r.text or "Nothing found" in r.text:
                    out.append(r.text)
                if "%J" in r.text or "Nothing found" in r.text:  # a definite answer: do not ask the other mirrors
                    break
        except httpx.HTTPError as e:
            log.warning("sesame %s failed: %s", base, e)
    if not out:
        raise ApiError(502, "upstream", "SIMBAD/Sesame недоступны")
    return out


def _simbad_section(text: str) -> list[str]:
    """Only the lines Sesame took from SIMBAD: NED/VizieR blocks may lack proper motion, never mix them in."""
    keep, on = [], False
    for line in text.splitlines():
        if line.startswith("#="):
            if on and keep:
                break
            on = "=Sc=Simbad" in line
        elif on:
            keep.append(line)
    return keep


def _parse(name: str, texts: list[str]) -> dict | None:
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    for text in texts:
        lines = _simbad_section(text)
        j = next((re.match(rf"%J\s+({number})\s+({number})(?:\s|$)", l) for l in lines if l.startswith("%J ")), None)
        if not j:
            if any(line.startswith("%J ") for line in lines):
                raise ApiError(502, "upstream", "Не удалось разобрать координаты в ответе Sesame")
            continue
        if not (0 <= float(j[1]) < 360 and -90 <= float(j[2]) <= 90):
            raise ApiError(502, "upstream", "Sesame вернул недопустимые координаты")
        pm = next((re.match(rf"%P\s+({number})\s+({number})(?:\s|$)", l) for l in lines if l.startswith("%P ")), None)
        x = next((re.match(rf"%X\s+({number})(?:\s|$)", l) for l in lines if l.startswith("%X ")), None)
        rv = next((re.match(rf"%V\s+v\s+({number})(?:\s|$)", l) for l in lines if l.startswith("%V v ")), None)
        if any(not math.isfinite(float(value)) for match in (pm, x, rv) if match for value in match.groups()):
            raise ApiError(502, "upstream", "Sesame вернул нечисловые параметры движения")
        return {"name": name, "ra": float(j[1]), "dec": float(j[2]),
                "pmra_masyr": float(pm[1]) if pm else None, "pmdec_masyr": float(pm[2]) if pm else None,
                "parallax_mas": float(x[1]) if x else None, "ref_epoch": None,
                "source": "SIMBAD (via CDS Sesame)", "resolver_version": 3, "ref_epoch_format": None,
                "ref_epoch_scale": None, "source_id": next((l.split(maxsplit=1)[1].strip() for l in lines if l.startswith("%I") and len(l.split(maxsplit=1)) == 2), None),
                "source_url": "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?" + quote(name, safe=""),
                "rv_kms": float(rv[1]) if rv else None, "errors": None, "coordinate_model": "fixed_catalogue",
                "assumptions": ["reference_epoch_unverified", "proper_motion_not_applied"],
                "limitations": "Sesame J2000 position; reference epoch and time scale are not verified."}
    return None


_team = {}
PUBLIC_IDS = {"barnard": "barnard", "epsilon_indi_b": "epsindib", "wise_0855": "wise0855",
              "rw_aur": "rwaur", "r_dra": "rdra"}


def _catalog():
    import sky_motion
    stamp = sky_motion.DEFAULT_CATALOG.stat().st_mtime_ns
    if _team.get("stamp") != stamp:
        _team.update(stamp=stamp, targets=sky_motion.load_targets())
    return _team["targets"]


def _preset_record(query: str) -> dict | None:
    from processing.targets import find_preset, load_targets
    norm = lambda value: "".join(c for c in value.lower() if c.isalnum())
    alias = next((p["id"] for p in PRESETS if norm(p["name"]) == norm(query)), query)
    match = find_preset(PUBLIC_IDS.get(alias, alias), load_targets())
    if match is None:
        return None
    t = _catalog()[match.name]
    row = t.row
    def number(key):
        return float(row[key]) if row.get(key, "").strip() else None
    return {"name": t.name, "id": match.id, "ra": number("ra_deg"), "dec": number("dec_deg"),
            "pmra_masyr": number("pmra_mas_yr"), "pmdec_masyr": number("pmdec_mas_yr"),
            "parallax_mas": number("parallax_mas"), "rv_kms": number("radial_velocity_km_s"),
            "ref_epoch": number("ref_epoch"), "ref_epoch_format": row["ref_epoch_format"],
            "ref_epoch_scale": row["ref_epoch_scale"].lower() or None,
            "source": row["catalog"], "source_id": row["source_id"], "source_url": row["source_url"],
            "canonical_id": row["canonical_id"], "component": row["component"],
            "coordinate_frame": row["coordinate_frame"], "ra_error_basis": row["ra_error_basis"],
            "coordinate_model": "geocenter_geometric", "assumptions": list(t.assumptions),
            "limitations": row["notes"], "errors": {k: number(k) for k in
                ("ra_error_mas", "dec_error_mas", "pmra_error_mas_yr", "pmdec_error_mas_yr",
                 "parallax_error_mas", "radial_velocity_error_km_s")}}


def resolve(name: str) -> dict:
    """Preset id or name -> team table; any other name -> SIMBAD (cached forever: memory + resolve.json).
    Result: position at ref_epoch, proper motion, parallax."""
    rec = _preset_record(name.strip())
    if rec:
        return rec
    k = name.strip().lower()
    with _lock:
        if not _mem and _FILE.exists():
            _mem.update(json.loads(_FILE.read_text(encoding="utf-8")))
        if k in _mem and _mem[k].get("resolver_version") == 3:
            return _mem[k]
    rec = _parse(name.strip(), _fetch_sesame(name.strip()))
    if rec is None:
        raise ApiError(404, "not_found", f"объект {name!r} не найден в SIMBAD")
    with _lock:
        _mem[k] = rec
        _FILE.write_text(json.dumps(_mem, ensure_ascii=False), encoding="utf-8")
    return rec


def propagate(o: dict, mjd: float):
    """MJD UTC -> common geocentric prediction, or space motion from a known reference epoch; otherwise fixed."""
    if o.get("coordinate_model") == "geocenter_geometric":
        import sky_motion
        pos = sky_motion.predict_position(_catalog()[o["name"]], Time(mjd, format="mjd", scale="utc"))
        return float(pos.geocentric.ra.deg), float(pos.geocentric.dec.deg)
    if o.get("ref_epoch") is None or o.get("pmra_masyr") is None or o.get("pmdec_masyr") is None:
        return o["ra"], o["dec"]
    kw = dict(ra=o["ra"] * u.deg, dec=o["dec"] * u.deg,
              pm_ra_cosdec=o["pmra_masyr"] * u.mas / u.yr, pm_dec=o["pmdec_masyr"] * u.mas / u.yr,
              obstime=Time(o["ref_epoch"], format=o.get("ref_epoch_format", "jyear"), scale=o.get("ref_epoch_scale") or "tt"))
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
    mjd = center_mjd()
    return {"center": propagate(o, mjd), "motion": o, "label": o["name"], "center_mjd": mjd}
