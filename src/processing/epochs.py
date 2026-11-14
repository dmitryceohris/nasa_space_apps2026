"""Группировка снимков во «временные эпохи» — кадры будущего видео."""
from __future__ import annotations

import numpy as np


def group_epochs(mjds, gap_days: float = 30.0, max_span_days: float = 20.0) -> list[list[int]]:
    """Индексы, сгруппированные по эпохам.

    1. Новая эпоха начинается, если между соседними снимками > gap_days
       (SPHEREx проходит каждую точку неба раз в ~6 месяцев).
    2. Если эпоха длиннее max_span_days (поля у полюсов эклиптики наблюдаются почти непрерывно),
       она режется на окна по max_span_days — иначе «видео» у R Dra было бы из одного кадра.
    """
    mjds = np.asarray(mjds, dtype=float)
    if mjds.size == 0:
        return []
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
