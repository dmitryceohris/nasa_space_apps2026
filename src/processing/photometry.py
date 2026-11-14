"""Фотометрия в заданной точке: спектр «на лету» и кривая блеска.

Измеряем на РОДНЫХ пикселях каждой экспозиции (без пересчёта на сетку — интерполяция
искажает поток). Длина волны точки берётся из спектрального WCS экспозиции: из-за
LVF-фильтра каждая экспозиция видит объект на своей λ — так и собирается спектр.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np
from astropy.coordinates import SkyCoord

from Exposure import Exposure

PHOTOMETRY_METHOD = {
    "name": "forced_aperture_photometry",
    "description": "Круглая апертура r=2 пикс. на родных пикселях каждой экспозиции (без перепроецирования), "
                   "фон — медиана кольца 4–7 пикс. после вычитания модели зодиакального света; "
                   "положение цели — каталожное с учётом собственного движения, а не измеренное.",
    "params": {"r_ap_px": 2.0, "r_in_px": 4.0, "r_out_px": 7.0},
    "quality": "ok=false: в апертуре есть плохие (в т.ч. насыщенные) пиксели, поток восстановлен по хорошим",
    "limitations": [
        "без поправки на апертуру: часть потока PSF за пределами r=2 пикс. теряется (зависит от λ)",
        "соседние источники в апертуре или кольце искажают поток",
        "яркие цели насыщены: такие точки помечены ok=false и не входят в плоские массивы",
    ],
    "not_implemented": ["PSF-фотометрия", "апертурная поправка по λ"],
}


def forced_photometry(exp: Exposure, coord: SkyCoord, r_ap: float = 2.0,
                      r_in: float = 4.0, r_out: float = 7.0) -> dict | None:
    """Апертурная фотометрия: круг r_ap пикселей минус медиана кольца (r_in..r_out).

    Возвращает поток в μJy. ok=False, если в апертуре есть плохие пиксели (например,
    насыщение у ярких звёзд) — такие точки показываем, но помечаем.
    """
    if not exp.contains(coord, margin_px=0):
        return None
    x0, y0 = exp.sky_to_pixel(coord)
    yy, xx = np.indices(exp.shape)
    r = np.hypot(xx - x0, yy - y0)
    ap = r <= r_ap
    ann = (r >= r_in) & (r <= r_out)
    bad = exp.bad_mask()
    img = exp.background_subtracted()

    ann_good = ann & ~bad
    if ap.sum() == 0 or ann_good.sum() < 5:
        return None
    sky = float(np.median(img[ann_good]))
    ap_good = ap & ~bad
    n_ap, n_good = int(ap.sum()), int(ap_good.sum())
    if n_good == 0:
        flux = var = np.nan
    else:
        # недостающие (плохие) пиксели апертуры заполняем средним по хорошим
        flux = float(np.sum(img[ap_good] - sky)) * n_ap / n_good
        var = float(np.sum(exp.variance[ap_good])) * (n_ap / n_good) ** 2
    k = Exposure.mjysr_to_ujy_per_pixel()
    wl, bw = exp.wavelength_at(x0, y0)
    flags = sorted({n for yx in zip(*np.nonzero(ap)) for n in exp.flag_names_at(yx[1], yx[0])})
    return {"wave_um": round(wl, 5), "bw_um": round(bw, 5),
            "flux_ujy": None if not np.isfinite(flux) else round(flux * k, 3),
            "err_ujy": None if not np.isfinite(var) else round(float(np.sqrt(max(var, 0))) * k, 3),
            "mjd": round(exp.mjd, 5), "det": exp.detector, "obs_id": exp.obs_id,
            "ok": bool(n_good == n_ap), "flags": flags}


def spectrum(exps: list[Exposure], coord_at: Callable[[float], SkyCoord]) -> list[dict]:
    """Точки спектра со всех экспозиций. coord_at(mjd) учитывает собственное движение цели."""
    pts = [p for e in exps if (p := forced_photometry(e, coord_at(e.mjd))) is not None]
    return sorted(pts, key=lambda p: p["wave_um"])


def light_curve(points: list[dict], wave_min: float, wave_max: float) -> list[dict]:
    """Поток по времени в узком окне λ: иначе переменность спутается с формой спектра (LVF)."""
    sel = [p for p in points if wave_min <= p["wave_um"] <= wave_max and p["flux_ujy"] is not None]
    return sorted(sel, key=lambda p: p["mjd"])
