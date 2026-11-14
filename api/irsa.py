"""IRSA access: which SPHEREx exposures cover a point, and FITS cutouts. Every request is logged."""
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

from .cache import CACHE
from .errors import ApiError

log = logging.getLogger("irsa")
IRSA = "https://irsa.ipac.caltech.edu"
SEM = threading.BoundedSemaphore(6)  # at most 6 simultaneous requests to IRSA for the whole server
warnings.filterwarnings("ignore")


def _send(method, url, retries=6, **kw):
    """Request with backoff (IRSA answers 503 'Server busy' under load); logs url, time, status."""
    code = 0
    for i in range(retries):
        t = time.time()
        r = None
        with SEM:
            try:
                r = httpx.request(method, url, timeout=300, follow_redirects=True, **kw)
                code = r.status_code
            except httpx.HTTPError as e:
                code = 0
                log.warning("%s %s failed: %s", method, url, e)
        log.info("%s %s -> %s %.1fs", method, url, code, time.time() - t)
        if r is not None and code == 200:
            return r
        time.sleep(5 * (i + 1))
    raise ApiError(502, "upstream", f"IRSA недоступен (последний код {code})")


# ---------- local IPAC tables: a swappable epoch source ----------
# Any IRSA "Spectral Images" .tbl export dropped into TABLES is used when it covers the point.
TABLES = Path(os.environ.get("SPHEREX_TABLES", Path(__file__).resolve().parent.parent / "tables"))
_parsed = {}  # path -> (mtime, rows)


def _table_rows():
    rows = []
    for p in sorted(TABLES.glob("*.tbl")):
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
    for p, (_, rows) in _parsed.items():
        if p.exists() and rows:
            pts = np.vstack([r["poly"] for r in rows])
            out.append(dict(file=p.name, rows=len(rows), ra=[float(pts[:, 0].min()), float(pts[:, 0].max())],
                            dec=[float(pts[:, 1].min()), float(pts[:, 1].max())]))
    return out


def obs_id(uri: str) -> str | None:
    m = re.search(r"level2_(\d{4}W\d+_\w+?_\d+_\d)D\d_spx", uri) or re.search(r"\?(\d{4}W\w+)/D\d$", uri)
    return m[1] if m else None


# ---------- exposures covering a point ----------
def _tap_exposures(ra, dec, det):
    f = CACHE / f"exp_{ra:.3f}_{dec:.3f}_d{det}.json"
    if f.exists():
        return json.loads(f.read_text())
    q = ("SELECT a.uri,p.time_bounds_lower,p.energy_bounds_lower,p.energy_bounds_upper "
         "FROM spherex.plane p JOIN spherex.artifact a ON a.planeid=p.planeid "
         f"WHERE CONTAINS(POINT('ICRS',{ra},{dec}),p.poly)=1 AND p.dataproducttype='image' "
         f"AND p.energy_bandpassname='SPHEREx-D{det}'")
    r = _send("POST", f"{IRSA}/TAP/sync", data=dict(QUERY=q, FORMAT="csv"))
    out = []
    for line in r.text.strip().splitlines()[1:]:
        uri, mjd, lo, hi = line.split(",")
        if "level2_" in uri:
            out.append(dict(uri=uri, mjd=float(mjd), lam=[float(lo) * 1e6, float(hi) * 1e6]))
    f.write_text(json.dumps(out))
    return out


def exposures(ra, dec, det):
    """Flat list of {uri, mjd, lam:[lo,hi] um} for one detector: local table if it covers the point, else TAP."""
    out = [dict(uri=r["uri"], mjd=r["mjd"], lam=r["lam"]) for r in _table_rows()
           if r["band"] == f"SPHEREx-D{det}" and _inside(r["poly"], ra, dec)]
    return out or _tap_exposures(ra, dec, det)


def select(exps, max_per_epoch, gap_days=30.0):
    """Split by time gaps > gap_days into epochs, keep max_per_epoch exposures per epoch, evenly spaced in time."""
    exps = sorted(exps, key=lambda e: e["mjd"])
    clusters, cur = [], []
    for e in exps:
        if cur and e["mjd"] - cur[-1]["mjd"] > gap_days:
            clusters.append(cur)
            cur = []
        cur.append(e)
    if cur:
        clusters.append(cur)
    out = []
    for c in clusters:
        if len(c) > max_per_epoch:
            idx = np.unique(np.linspace(0, len(c) - 1, max_per_epoch).round().astype(int))
            c = [c[i] for i in idx]
        out.append(c)
    return out


# ---------- cutouts ----------
def _fits_url(uri: str) -> str:
    """TAP gives a relative IBE path; tables give a datalink URL whose '#this' row is the FITS."""
    if not uri.startswith("http"):
        return f"{IRSA}/{uri}"
    p = CACHE / f"dl_{hashlib.md5(uri.encode()).hexdigest()}.txt"
    if not p.exists():
        t = parse_single_table(io.BytesIO(_send("GET", uri).content)).to_table()
        p.write_text(str(next(r["col_1"] for r in t if r["col_4"] == "#this")))  # col_1=access_url, col_4=semantics
    return p.read_text()


def cutout(uri: str, ra: float, dec: float, size_deg: float) -> Path:
    h = hashlib.md5(f"{uri}{ra:.4f}{dec:.4f}{size_deg:.3f}".encode()).hexdigest()
    p = CACHE / f"cut_{h}.fits"
    if not p.exists():
        r = _send("GET", _fits_url(uri), params=dict(center=f"{ra},{dec}", size=f"{size_deg}deg", gzip="false"))
        p.write_bytes(r.content)
    return p
