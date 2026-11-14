"""Getting data: targets, archive search, downloads, local files, epochs."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from Exposure import Exposure, download_rows, load_exposure, query_spherex
from processing.epochs import GAP_DAYS, MAX_SPAN_DAYS, group_epochs, select_per_epoch
from processing.targets import Target, find_preset, load_targets, resolve_simbad

FIXED_EPOCH = Time(2016.0, format="jyear")  # reference epoch of a fixed (no motion) point; irrelevant for it


def target(where) -> Target:
    """Anything that names a place on the sky -> Target.

    'barnard', "Barnard's Star", 'wise0855' (the five challenge targets, offline) or any SIMBAD name (network);
    (ra_deg, dec_deg); SkyCoord; a Target (returned as is). Coordinates give a fixed point without motion.
    t.at(mjd) is the position on a date (proper motion + space motion; annual parallax is added by track())."""
    if isinstance(where, Target):
        return where
    if isinstance(where, SkyCoord):
        c = where.icrs
        return Target("coords", f"{c.ra.deg:.5f} {c.dec.deg:+.5f}", float(c.ra.deg), float(c.dec.deg),
                      0.0, 0.0, 0.0, 0.0, FIXED_EPOCH, source="coordinates")
    if isinstance(where, (tuple, list, np.ndarray)) and len(where) == 2:
        return target(SkyCoord(float(where[0]) * u.deg, float(where[1]) * u.deg))
    if isinstance(where, str):
        t = find_preset(where, load_targets()) or resolve_simbad(where)
        if t is None:
            raise ValueError(f"объект не найден: {where!r}")
        return t
    raise TypeError("цель: имя, (ra, dec) в градусах, SkyCoord или Target")


def _now_mjd() -> float:
    return float(Time.now().mjd)


def search(where, det: int | None = None) -> list[dict]:
    """Which SPHEREx exposures cover the target, without downloading. Needs network.
    -> [{obs_id, det, mjd, date, wave_min_um, wave_max_um, url}], sorted by time."""
    t = target(where)
    tbl = query_spherex(t.at(_now_mjd()), detector=det)
    out = []
    for r in tbl:
        out.append({"obs_id": str(r["obs_id"]), "det": int(r["detector"]), "mjd": float(r["mjd"]),
                    "date": Time(float(r["mjd"]), format="mjd").to_value("iso", subfmt="date"),
                    "wave_min_um": float(r["em_min"]) * 1e6, "wave_max_um": float(r["em_max"]) * 1e6,
                    "url": str(r["access_url"])})
    return out


def fetch(where, det: int = 3, fov: float = 3.0, per_epoch: int = 6, progress=None) -> list[Exposure]:
    """Search IRSA and download cutouts around the target: at most `per_epoch` exposures per epoch
    (same choice as the web site). Needs network; cutouts are cached by astropy.
    progress(done, total) is called during the download."""
    t = target(where)
    center = t.at(_now_mjd())
    tbl = query_spherex(center, detector=det)
    if len(tbl) == 0:
        return []
    sel = [i for g in select_per_epoch(np.asarray(tbl["mjd"], float), per_epoch) for i in g]
    size_deg = fov / 60 * 1.45 + 0.01  # the north-up grid is rotated relative to the detector: take a margin
    exps = download_rows([tbl[i] for i in sel], center, size_deg, progress=progress)
    if len(exps) < len(sel):  # the loader skips failed downloads (IRSA answers 503 under load)
        import warnings

        warnings.warn(f"spxhack.fetch: скачано {len(exps)} из {len(sel)} снимков; повторный вызов докачает "
                      "остальные (готовые берутся из кэша)", stacklevel=2)
    return exps


def load(paths, **kw) -> Exposure | list[Exposure]:
    """A FITS path or URL -> Exposure; a list of them -> list. kw: mjd_hint, obs_id_hint, cache, timeout."""
    if isinstance(paths, (str, Path)):
        return load_exposure(str(paths), **kw)
    return [load_exposure(str(p), **kw) for p in paths]


def epochs(mjds, gap: float = GAP_DAYS, max_span: float = MAX_SPAN_DAYS) -> list[list[int]]:
    """Indices grouped into epochs: a gap > `gap` days starts a new one, none longer than `max_span` days."""
    seq = list(mjds)
    if seq and isinstance(seq[0], Exposure):
        seq = [e.mjd for e in seq]
    return group_epochs(seq, gap, max_span)
