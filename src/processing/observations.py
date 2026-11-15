"""Canonical product identity and provenance, independent of a particular selection."""
from __future__ import annotations

import hashlib
import re
from urllib.parse import urlsplit, urlunsplit


def product_record(url, *, obs_id="", detector=None, selection=None):
    parts = urlsplit(url)
    product_url = urlunsplit((parts.scheme, parts.netloc, parts.path, "", ""))
    match = re.search(r"/spherex/(dr1|qr[123])/level2/[^/]+/([^/]+)/([1-6])/", parts.path)
    release, version, det = match.groups() if match else (None,None,detector)
    if detector is not None and det is not None and int(det) != int(detector):
        raise ValueError("product_detector_mismatch")
    if not obs_id:
        name = re.search(r"level2_(.+)D[1-6]_spx_",parts.path)
        obs_id = name.group(1) if name else ""
    return {"product_id": hashlib.sha256(product_url.encode()).hexdigest()[:24],
            "product_url": product_url, "obs_id": obs_id, "detector": int(det) if det else None,
            "release": release, "processing_version": version,
            "selections": [selection] if selection else []}


def merge_registry(rows):
    """Selections refer to products; files/cutouts keep their own hashes and timestamps."""
    products = {}
    for row in rows:
        record = product_record(row["product_url"], obs_id=row.get("obs_id", ""),
                                detector=row.get("detector"), selection=row.get("selection"))
        key = record["product_id"]
        entry = products.setdefault(key, record | {"files": []})
        if (entry["obs_id"],entry["detector"]) != (record["obs_id"],record["detector"]):
            raise ValueError("conflicting_product_identity")
        entry["selections"] = sorted(set(entry["selections"]+record["selections"]))
        if row.get("sha256"):
            file = {k:row[k] for k in ("sha256","cutout_url","mjd","time_source") if k in row}
            if file not in entry["files"]:
                entry["files"].append(file)
    return sorted(products.values(),key=lambda r:r["product_url"])
