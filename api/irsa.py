"""IRSA access: which SPHEREx exposures cover a point, and FITS cutouts. Every request is logged."""
from contextlib import contextmanager
import hashlib
import io
import json
import logging
import os
import re
import threading
import time
import warnings
from pathlib import Path

import httpx
import numpy as np
from astropy.io import ascii
from astropy.io.votable import parse_single_table

from . import cache
from .cache import CACHE
from .errors import ApiError

log = logging.getLogger("irsa")
IRSA = "https://irsa.ipac.caltech.edu"
SEM = threading.BoundedSemaphore(6)  # at most 6 simultaneous requests to IRSA for the whole server
_inflight_lock = threading.Lock()
_inflight: dict[str, threading.Lock] = {}
warnings.filterwarnings("ignore")


@contextmanager
def _key_lock(k: str):
    with _inflight_lock:
        lk = _inflight.setdefault(k, threading.Lock())
    with lk:
        yield


def _irsa_base() -> str:
    return os.environ.get("IRSA_BASE_URL", IRSA).rstrip("/")


def _send(method, url, retries=None, **kw):
    """Request with exponential backoff (IRSA answers 503 'Server busy' under load); logs url, time, status."""
    if retries is None:
        retries = int(os.environ.get("IRSA_RETRIES", "6"))
    backoff = float(os.environ.get("IRSA_BACKOFF", "1.0"))
    timeout = float(os.environ.get("IRSA_TIMEOUT", "300.0"))
    code = 0
    for i in range(max(1, retries)):
        t = time.time()
        r = None
        with SEM:
            try:
                r = httpx.request(method, url, timeout=timeout, follow_redirects=True, **kw)
                code = r.status_code
            except httpx.HTTPError as e:
                code = 0
                log.warning("%s %s failed: %s", method, url, e)
        log.info("%s %s -> %s %.1fs", method, url, code, time.time() - t)
        if r is not None and code == 200:
            return r
        if i + 1 < retries:
            time.sleep(backoff * (2 ** i))
    raise ApiError(502, "upstream", f"IRSA недоступен (последний код {code})")


# ---------- local IPAC tables: a swappable epoch source ----------
# Any IRSA "Spectral Images" .tbl export dropped into TABLES is used when it covers the point.
_DEFAULT_TABLES = Path(__file__).resolve().parent.parent / "tables"
_parsed = {}  # path -> (mtime, rows)


def _tables_dir() -> Path:
    return Path(os.environ.get("SPHEREX_TABLES", _DEFAULT_TABLES))


def _table_rows():
    rows = []
    tdir = _tables_dir()
    if not tdir.exists():
        return rows
    for p in sorted(tdir.glob("*.tbl")):
        m = p.stat().st_mtime
        if p not in _parsed or _parsed[p][0] != m:
            try:
                t = ascii.read(p, format="ipac")
                rows_p = [dict(uri=str(r["access_url"]), band=str(r["energy_bandpassname"]), mjd=float(r["t_min"]),
                               lam=[float(r["em_min"]) * 1e6, float(r["em_max"]) * 1e6],
                               poly=np.array(str(r["s_region"]).split()[2:], float).reshape(-1, 2)) for r in t]
            except Exception as e:  # not an exposure table (e.g. a column-description table): skip it
                log.warning("skipping %s: %s", p.name, e)
                rows_p = []
            _parsed[p] = (m, rows_p)
        rows += _parsed[p][1]
    return rows


def _inside(poly, ra, dec):
    """Ray casting in the tangent plane around the point. ponytail: flat approximation, fine for ~3 deg
    footprints away from the poles; use a spherical polygon test if polar targets matter."""
    x = ((poly[:, 0] - ra + 180) % 360 - 180) * np.cos(np.radians(dec))
    y = poly[:, 1] - dec
    ins = False
    for i in range(len(x) - 1):
        if (y[i] > 0) != (y[i + 1] > 0) and x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]) > 0:
            ins = not ins
    return ins


def sources() -> list:
    """Loaded tables with row count and sky bounding box."""
    _table_rows()
    out = []
    tdir = _tables_dir()
    for p, (_, rows) in _parsed.items():
        if p.parent == tdir and p.exists() and rows:
            pts = np.vstack([r["poly"] for r in rows])
            out.append(dict(file=p.name, rows=len(rows), ra=[float(pts[:, 0].min()), float(pts[:, 0].max())],
                            dec=[float(pts[:, 1].min()), float(pts[:, 1].max())]))
    return out


def obs_id(uri: str) -> str | None:
    m = re.search(r"level2_(\d{4}W\d+_\w+?_\d+_\d)D\d_spx", uri) or re.search(r"\?(\d{4}W\w+)/D\d$", uri)
    return m[1] if m else None


# ---------- exposures covering a point ----------
def _tap_exposures(ra, dec, det):
    name = f"exp_{ra:.3f}_{dec:.3f}_d{det}.json"
    with _key_lock(name):
        f = cache.CACHE / name
        ttl = float(os.environ.get("IRSA_EXP_TTL", "86400"))
        stale_data = None
        if f.exists():
            try:
                stale_data = json.loads(f.read_text(encoding="utf-8"))
                if time.time() - f.stat().st_mtime < ttl:
                    return stale_data
            except Exception:
                stale_data = None
        q = ("SELECT a.uri,p.time_bounds_lower,p.energy_bounds_lower,p.energy_bounds_upper "
             "FROM spherex.plane p JOIN spherex.artifact a ON a.planeid=p.planeid "
             f"WHERE CONTAINS(POINT('ICRS',{ra},{dec}),p.poly)=1 AND p.dataproducttype='image' "
             f"AND p.energy_bandpassname='SPHEREx-D{det}'")
        try:
            r = _send("POST", f"{_irsa_base()}/TAP/sync", data=dict(QUERY=q, FORMAT="csv"))
        except ApiError:
            if stale_data is not None:
                log.warning("TAP unavailable; using expired exposure list %s", name)
                return stale_data
            raise
        out = []
        for line in r.text.strip().splitlines()[1:]:
            uri, mjd, lo, hi = line.split(",")
            if "level2_" in uri:
                out.append(dict(uri=uri, mjd=float(mjd), lam=[float(lo) * 1e6, float(hi) * 1e6]))
        cache.atomic_write_text(f, json.dumps(out))
        return out


def exposures(ra, dec, det):
    """Flat list of {uri, mjd, lam:[lo,hi] um} for one detector: local table if it covers the point, else TAP."""
    out = [dict(uri=r["uri"], mjd=r["mjd"], lam=r["lam"]) for r in _table_rows()
           if r["band"] == f"SPHEREx-D{det}" and _inside(r["poly"], ra, dec)]
    return out or _tap_exposures(ra, dec, det)


# ---------- cutouts ----------
def _fits_url(uri: str) -> str:
    """TAP gives a relative IBE path; tables give a datalink URL whose '#this' row is the FITS."""
    if not uri.startswith("http"):
        return f"{_irsa_base()}/{uri.lstrip('/')}"
    name = f"dl_{hashlib.md5(uri.encode()).hexdigest()}.txt"
    with _key_lock(name):
        p = cache.CACHE / name
        if not p.exists():
            t = parse_single_table(io.BytesIO(_send("GET", uri).content)).to_table()
            url = str(next(r["col_1"] for r in t if r["col_4"] == "#this"))  # col_1=access_url, col_4=semantics
            cache.atomic_write_text(p, url)
        return p.read_text(encoding="utf-8")


def cutout(uri: str, ra: float, dec: float, size_deg: float) -> Path:
    h = hashlib.md5(f"{uri}{ra:.4f}{dec:.4f}{size_deg:.3f}".encode()).hexdigest()
    name = f"cut_{h}.fits"
    with _key_lock(name):
        p = cache.CACHE / name
        if not p.exists():
            r = _send("GET", _fits_url(uri), params=dict(center=f"{ra},{dec}", size=f"{size_deg}deg", gzip="false"))
            cache.atomic_write_bytes(p, r.content)
        return p
