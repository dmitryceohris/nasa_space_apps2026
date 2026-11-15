"""IRSA access: which SPHEREx exposures cover a point, and FITS cutouts. Every request is logged."""
import hashlib
import csv
import io
import json
import logging
import os
import re
import threading
import tempfile
import time
import warnings
from pathlib import Path
from urllib.parse import urlsplit

import httpx
import numpy as np
from astropy.io import ascii, fits
from astropy.io.fits.verify import VerifyError
from astropy.io.votable import parse_single_table

from Exposure.irsa import archive_collections

from .cache import CACHE
from .errors import ApiError
from runtime_storage import current_workspace

log = logging.getLogger("irsa")
IRSA = "https://irsa.ipac.caltech.edu"
SEM = threading.BoundedSemaphore(6)  # at most 6 simultaneous requests to IRSA for the whole server
_tap_probe_lock = threading.Lock()
_tap_retry_at = 0.0
SIA_LIMITATION = "Поиск выполнен через IRSA SIA2 (текущие коллекции): новые наблюдения могут появляться примерно на сутки позже TAP"
warnings.filterwarnings("ignore")


def _send(method, url, retries=3, **kw):
    """Request with backoff (IRSA answers 503 'Server busy' under load); logs url, time, status."""
    code = 0
    timeout = kw.pop("timeout", 45)
    for i in range(retries):
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
        if 400 <= code < 500 and code not in (408, 429):  # the archive answered: retrying will not help
            raise ApiError(502, "upstream", f"IRSA отклонил запрос (код {code})")
        if i + 1 < retries:
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
class ExposureRows(list):
    """Exposure records with provenance, including the empty-result case."""

    def __init__(self, rows=(), *, source="tap", warnings=()):
        super().__init__(rows)
        self.source = source
        self.warnings = list(warnings)


def _sia_exposures(ra, dec, det):
    collections = archive_collections()
    response = _send("GET", f"{IRSA}/SIA", params={
        "POS": f"CIRCLE {ra} {dec} 0", "COLLECTION": list(collections),
        "DPTYPE": "image", "MAXREC": "10000",
    })
    try:
        from astropy.io.votable import parse

        document = parse(io.BytesIO(response.content))
        for resource in document.resources:
            for info in resource.infos:
                if info.name == "QUERY_STATUS" and info.value != "OK":
                    raise ValueError(f"SIA QUERY_STATUS={info.value}: {info.content}")
        table = document.get_first_table().to_table(use_names_over_ids=True)
        rows, seen = [], set()
        # DR1 takes priority over duplicate QR3 entries, matching the Python library.
        rank = {collection: i for i, collection in enumerate(collections)}
        for row in sorted(table, key=lambda r: rank.get(str(r["obs_collection"]), 99)):
            if str(row["energy_bandpassname"]) != f"SPHEREx-D{det}":
                continue
            identity = str(row["obs_id"])
            if identity in seen:
                continue
            seen.add(identity)
            rows.append({"uri": str(row["access_url"]),
                         "mjd": (float(row["t_min"]) + float(row["t_max"])) / 2,
                         "lam": [float(row["em_min"]) * 1e6, float(row["em_max"]) * 1e6]})
        if any(not np.isfinite([r["mjd"], *r["lam"]]).all() for r in rows):
            raise ValueError("SIA returned non-finite metadata")
        return ExposureRows(sorted(rows, key=lambda r: r["mjd"]), source="sia2", warnings=[SIA_LIMITATION])
    except (ValueError, KeyError, IndexError, TypeError) as error:
        log.warning("invalid SIA response: %s", error)
        raise ApiError(502, "upstream", "IRSA SIA2 вернул непригодный или неполный ответ") from error


def _tap_exposures(ra, dec, det):
    global _tap_retry_at
    release = hashlib.sha256(",".join(archive_collections()).encode()).hexdigest()[:12]
    f = CACHE / f"exp_{ra:.3f}_{dec:.3f}_d{det}_{release}.json"
    try:
        saved = json.loads(f.read_text())
        if isinstance(saved, dict):
            saved = ExposureRows(saved["rows"], source=saved["source"], warnings=saved["warnings"])
        elif isinstance(saved, list):
            saved = ExposureRows(saved)
        else:
            raise ValueError("invalid exposure cache")
        ttl = float(os.environ.get("EXPOSURE_CACHE_TTL", "3600")) if saved else float(os.environ.get("EMPTY_CACHE_TTL", "60"))
        if time.time() - f.stat().st_mtime < ttl:
            return saved
    except (OSError, ValueError, KeyError, TypeError):
        pass
    q = ("SELECT a.uri,p.time_bounds_lower,p.energy_bounds_lower,p.energy_bounds_upper "
         "FROM spherex.plane p JOIN spherex.artifact a ON a.planeid=p.planeid "
         f"WHERE CONTAINS(POINT('ICRS',{ra},{dec}),p.poly)=1 AND p.dataproducttype='image' "
         f"AND p.energy_bandpassname='SPHEREx-D{det}'")
    out = None
    # One bounded probe for concurrent jobs; a broken TAP backend must not delay every detector.
    with _tap_probe_lock:
        if time.monotonic() >= _tap_retry_at:
            try:
                r = _send("GET", f"{IRSA}/TAP/sync", params=dict(QUERY=q, FORMAT="csv"),
                          retries=1, timeout=float(os.environ.get("TAP_PROBE_TIMEOUT", "12")))
                records = list(csv.reader(io.StringIO(r.text)))
                if not records or len(records[0]) != 4:
                    raise ValueError("invalid TAP CSV response")
                out = ExposureRows(dict(uri=uri, mjd=float(mjd), lam=[float(lo) * 1e6, float(hi) * 1e6])
                                   for uri, mjd, lo, hi in records[1:] if "level2_" in uri)
                _tap_retry_at = 0.0
            except (ApiError, ValueError) as error:
                log.warning("TAP probe failed; selecting official SIA2: %s", error)
                _tap_retry_at = time.monotonic() + float(os.environ.get("TAP_RETRY_SECONDS", "300"))
    if out is None:
        out = _sia_exposures(ra, dec, det)
    _atomic_write(f, json.dumps({"rows": out, "source": out.source, "warnings": out.warnings}, allow_nan=False).encode())
    return out


def selection_token():
    mode = os.environ.get("SPHEREX_SOURCE", "archive")
    if mode == "archive":
        return mode
    if mode != "tables":
        raise ValueError("SPHEREX_SOURCE must be archive or tables")
    return mode + str([(p.name, p.stat().st_mtime_ns) for p in sorted(TABLES.glob("*.tbl"))])


def exposures(ra, dec, det):
    """Live archive by default. SPHEREX_SOURCE=tables explicitly selects a limited local sample."""
    if selection_token() == "archive":
        return _tap_exposures(ra, dec, det)
    return [dict(uri=r["uri"], mjd=r["mjd"], lam=r["lam"]) for r in _table_rows()
            if r["band"] == f"SPHEREx-D{det}" and _inside(r["poly"], ra, dec)]


def _atomic_write(path, payload, workspace=None):
    workspace = workspace or current_workspace()
    if workspace is not None:
        workspace.atomic_write(path, [payload])
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(payload)
            stream.flush()
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


# ---------- cutouts ----------
def _fits_url(uri: str, workspace=None) -> str:
    """TAP gives a relative IBE path; tables give a datalink URL whose '#this' row is the FITS."""
    if not uri.startswith("http"):
        return f"{IRSA}/{uri}"
    if urlsplit(uri).path.lower().endswith(".fits"):
        return uri  # SIA2 supplies direct FITS URLs; only datalink URLs need XML lookup.
    p = CACHE / f"dl_{hashlib.md5(uri.encode()).hexdigest()}.txt"
    if not p.exists():
        t = parse_single_table(io.BytesIO(_send("GET", uri).content)).to_table()
        _atomic_write(p, str(next(r["col_1"] for r in t if r["col_4"] == "#this")).encode(), workspace)
    return p.read_text()


def cutout(uri: str, ra: float, dec: float, size_deg: float, *, workspace=None) -> Path:
    h = hashlib.md5(f"{uri}{ra:.4f}{dec:.4f}{size_deg:.3f}".encode()).hexdigest()
    p = (workspace.path if workspace is not None else CACHE) / f"cut_{h}.fits"
    try:
        with fits.open(p) as hdul:
            hdul.verify("exception")
            if hdul["IMAGE"].data is None or hdul["IMAGE"].data.ndim != 2:
                raise ValueError("missing image")
        return p
    except (OSError, KeyError, ValueError, TypeError, VerifyError):
        pass
    r = _send("GET", _fits_url(uri, workspace), params=dict(center=f"{ra},{dec}", size=f"{size_deg}deg", gzip="false"))
    try:
        with fits.open(io.BytesIO(r.content)) as hdul:
            hdul.verify("exception")
            if hdul["IMAGE"].data is None or hdul["IMAGE"].data.ndim != 2:
                raise ValueError("missing image")
    except (OSError, KeyError, ValueError, TypeError, VerifyError) as error:
        raise ApiError(502, "upstream", "IRSA вернул повреждённый FITS") from error
    _atomic_write(p, r.content, workspace)
    return p
