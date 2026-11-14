"""Группировка снимков во «временные эпохи» — кадры будущего видео."""
from __future__ import annotations

import numpy as np

GAP_DAYS = 30.0       # разрыв между снимками, после которого начинается новая эпоха
MAX_SPAN_DAYS = 20.0  # самая длинная эпоха; длиннее — режем на окна


def epoch_limits(gap_days: float = GAP_DAYS, max_span_days: float = MAX_SPAN_DAYS) -> dict:
    """Ограничения эпох для ответа API: пользователь видит, какие снимки сложены вместе."""
    return {"gap_days": gap_days, "max_span_days": max_span_days,
            "description": "снимки с разрывом больше gap_days — разные эпохи; эпоха не длиннее max_span_days "
                           "(иначе быстрый объект размажется внутри одного кадра)"}


def smear_px(pm_arcsec_yr: float, span_days: float, pix_arcsec: float) -> float:
    """Насколько объект с собственным движением pm сдвигается за эпоху длиной span_days, в пикселях."""
    return abs(pm_arcsec_yr) * span_days / 365.25 / pix_arcsec


def group_epochs(mjds, gap_days: float = GAP_DAYS, max_span_days: float = MAX_SPAN_DAYS) -> list[list[int]]:
    """Индексы, сгруппированные по эпохам.

    1. Новая эпоха начинается, если между соседними снимками > gap_days
       (SPHEREx проходит каждую точку неба раз в ~6 месяцев).
    2. Если эпоха длиннее max_span_days (поля у полюсов эклиптики наблюдаются почти непрерывно),
       она режется на окна по max_span_days — иначе «видео» у R Dra было бы из одного кадра.
    Гарантии: каждый индекс ровно в одной эпохе, эпохи идут по времени, mjd_max − mjd_min ≤ max_span_days.
    """
    if not (gap_days > 0 and max_span_days >= 0):
        raise ValueError(f"нужно gap_days > 0 и max_span_days >= 0, получено {gap_days}, {max_span_days}")
    mjds = np.asarray(mjds, dtype=float)
    if mjds.size == 0:
        return []
    if not np.isfinite(mjds).all():
        raise ValueError("MJD снимков должны быть конечными числами")
    order = np.argsort(mjds)
    groups, cur = [], [int(order[0])]
    for prev, i in zip(order[:-1], order[1:]):
        if mjds[i] - mjds[prev] > gap_days:
            groups.append(cur)
            cur = []
        cur.append(int(i))
    groups.append(cur)

    out = []
    for g in groups:
        start = mjds[g[0]]
        chunk = []
        for i in g:
            if chunk and mjds[i] - start > max_span_days:
                out.append(chunk)
                chunk, start = [], mjds[i]
            chunk.append(i)
        out.append(chunk)
    return out


def pick_evenly(indices: list[int], n: int) -> list[int]:
    """n индексов, равномерно распределённых по списку (список уже отсортирован по времени)."""
    if n <= 0 or len(indices) <= n:
        return list(indices)
    pos = np.linspace(0, len(indices) - 1, n).round().astype(int)
    return [indices[p] for p in sorted(set(pos))]


def select_per_epoch(mjds, max_per_epoch: int, **kw) -> list[list[int]]:
    """Эпохи, в каждой не больше max_per_epoch снимков — чтобы не качать все ~60 вырезок."""
    return [pick_evenly(g, max_per_epoch) for g in group_epochs(mjds, **kw)]
