"""Поиск снимков SPHEREx в архиве NASA/IPAC IRSA и загрузка вырезок в Exposure.

Данные запрашиваются «вживую», поэтому новые релизы SPHEREx появляются в инструменте
автоматически (пункт челленджа «New SPHEREx data is arriving fast»).
"""
from __future__ import annotations

import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from .exposure import Exposure
from .loader import load_exposure

# Порядок важен: при дублях одной экспозиции в разных релизах берём первый (самый новый).
DEFAULT_COLLECTIONS = ("spherex_dr1", "spherex_qr3")


def query_spherex(coord: SkyCoord, *, detector: int | None = None,
                  collections=DEFAULT_COLLECTIONS):
    """Все спектральные снимки SPHEREx, покрывающие точку `coord`.

    Возвращает astropy Table со столбцами SIA (access_url, obs_id, t_min, t_max,
    energy_bandpassname, ...) + добавленными detector и mjd, отсортированную по времени.
    """
    from astroquery.ipac.irsa import Irsa

    tbl = Irsa.query_sia(pos=(coord, 1 * u.arcsec), collection=list(collections))
    if len(tbl) == 0:
        return tbl

    band = np.array([str(b) for b in tbl["energy_bandpassname"]])
    det = np.array([int(b.rsplit("D", 1)[-1]) if b.startswith("SPHEREx-D") else -1 for b in band])
    tbl["detector"] = det
    tbl["mjd"] = 0.5 * (np.asarray(tbl["t_min"], float) + np.asarray(tbl["t_max"], float))
    tbl = tbl[det > 0]
    if detector is not None:
        tbl = tbl[tbl["detector"] == detector]

    # Убираем дубли (одна и та же экспозиция в DR1 и QR3)
    rank = {c: i for i, c in enumerate(collections)}
    order = np.argsort([rank.get(str(c), 99) for c in tbl["obs_collection"]], kind="stable")
    tbl = tbl[order]
    seen, keep = set(), []
    for i, row in enumerate(tbl):
        key = (str(row["obs_id"]), int(row["detector"]))
        if key not in seen:
            seen.add(key)
            keep.append(i)
    tbl = tbl[keep]
    tbl.sort("mjd")
    return tbl


def cutout_url(access_url: str, coord: SkyCoord, size_deg: float) -> str:
    """URL серверной вырезки IRSA: те же HDU, но только квадрат size_deg вокруг coord."""
    return f"{access_url}?center={coord.ra.deg:.6f},{coord.dec.deg:.6f}d&size={size_deg:.5f}"


def download_rows(rows, coord: SkyCoord, size_deg: float, *, workers: int = 6,
                  progress=None) -> list[Exposure]:
    """Скачивает вырезки для выбранных строк таблицы query_spherex (параллельно).

    progress(done, total) вызывается после каждого файла. Сбойные загрузки пропускаются.
    """
    rows = list(rows)
    total, done = len(rows), 0

    def _one(row):
        try:
            return load_exposure(cutout_url(str(row["access_url"]), coord, size_deg),
                                 mjd_hint=float(row["mjd"]), obs_id_hint=str(row["obs_id"]))
        except Exception as e:  # noqa: BLE001 — сеть IRSA иногда рвётся, один кадр не важен
            warnings.warn(f"пропущен {row['obs_id']}: {e}")
            return None

    exps = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        # as_completed: прогресс растёт по мере готовности, а не в порядке очереди
        for fut in as_completed([pool.submit(_one, r) for r in rows]):
            e = fut.result()
            done += 1
            if progress:
                progress(done, total)
            if e is not None:
                exps.append(e)
    return sorted(exps, key=lambda e: e.mjd)


def fetch_exposures(coord: SkyCoord, size_deg: float = 0.1, *, detector: int | None = None,
                    max_files: int | None = None, workers: int = 6,
                    collections=DEFAULT_COLLECTIONS) -> list[Exposure]:
    """Ищет снимки и скачивает ВЫРЕЗКИ (не полные 70-МБ файлы). Сбойные загрузки пропускает."""
    tbl = query_spherex(coord, detector=detector, collections=collections)
    if max_files is not None:
        tbl = tbl[:max_files]
    return download_rows(tbl, coord, size_deg, workers=workers)
