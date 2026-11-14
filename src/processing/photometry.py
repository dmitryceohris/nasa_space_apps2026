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
    "quality": "ok=false: плохие/загрязнённые пиксели в апертуре, неполная апертура/фоновое кольцо или нечисловой результат",
    "limitations": [
        "без поправки на апертуру: часть потока PSF за пределами r=2 пикс. теряется (зависит от λ)",
        "соседние источники в апертуре или кольце искажают поток",
        "яркие цели насыщены: такие точки помечены ok=false и не входят в плоские массивы",
        "ошибка включает дисперсию пикселей апертуры (VARIANCE) и ошибку медианы фонового кольца; "
        "не включены: апертурная поправка, калибровка потока, ошибка положения цели (апертура ставится по "
        "каталожному прогнозу, а не по измеренному центру), ширина и форма спектрального канала",
    ],
    "error_components": ["aperture_pixel_variance", "background_median"],
    "error_not_included": ["aperture_correction", "flux_calibration", "position_error", "bandpass"],
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
    from .frames import lone_spike_mask

    img, bad = exp.background_subtracted(), exp.bad_mask()
    bad = bad | lone_spike_mask(exp, img, bad)
    contaminated = exp.bad_mask(include_contamination=True) & ~bad

    ann_good = ann & ~bad & ~contaminated
    if ap.sum() == 0 or ann_good.sum() < 5:
        return None
    sky_vals = img[ann_good]
    sky = float(np.median(sky_vals))
    # error of the median of the ring: (pi/2) * sigma^2 / N, sigma from the ring itself (MAD) so that faint
    # sources and structure in the ring count too
    sky_sig = 1.4826 * float(np.median(np.abs(sky_vals - sky)))
    var_sky_median = (np.pi / 2) * sky_sig ** 2 / ann_good.sum()
    ap_good = ap & ~bad
    n_ap, n_good = int(ap.sum()), int(ap_good.sum())
    if n_good == 0:
        flux = var = np.nan
    else:
        # недостающие (плохие) пиксели апертуры заполняем средним по хорошим
        flux = float(np.sum(img[ap_good] - sky)) * n_ap / n_good
        var = float(np.sum(exp.variance[ap_good])) * (n_ap / n_good) ** 2 + n_ap ** 2 * var_sky_median
    k = Exposure.mjysr_to_ujy_per_pixel()
    wl, bw = exp.wavelength_at(x0, y0)
    flags = sorted({n for yx in zip(*np.nonzero(ap)) for n in exp.flag_names_at(yx[1], yx[0])})
    aperture_complete = exp.contains(coord, margin_px=r_ap)
    annulus_complete = exp.contains(coord, margin_px=r_out)
    reasons = []
    if not aperture_complete:
        reasons.append("aperture_outside_cutout")
    if not annulus_complete:
        reasons.append("background_annulus_outside_cutout")
    if ann_good.sum() < 0.5 * ann.sum():
        reasons.append("background_annulus_mostly_masked")
    if n_good != n_ap:
        reasons.append("masked_aperture_pixels")
    if np.any(ap & contaminated):
        reasons.append("contaminated_aperture")
    if not np.isfinite(flux) or not np.isfinite(var):
        reasons.append("nonfinite_measurement")
    if not np.isfinite(wl) or not np.isfinite(bw):
        reasons.append("invalid_wavelength")
    return {"wave_um": round(wl, 5) if np.isfinite(wl) else None,
            "bw_um": round(bw, 5) if np.isfinite(bw) else None,
            "flux_ujy": None if not np.isfinite(flux) else round(flux * k, 3),
            "err_ujy": None if not np.isfinite(var) else round(float(np.sqrt(max(var, 0))) * k, 3),
            "mjd": round(exp.mjd, 5), "det": exp.detector, "obs_id": exp.obs_id,
            "ok": not reasons, "flags": flags, "quality_reasons": reasons,
            "coverage": {"aperture_complete": bool(aperture_complete),
                         "background_annulus_complete": bool(annulus_complete),
                         "aperture_good_pixels": n_good, "aperture_pixels": n_ap,
                         "background_good_pixels": int(ann_good.sum()),
                         "background_pixels": int(ann.sum())},
            "wavelength_source": "spectral_wcs" if exp.wave_wcs is not None else "approximate_detector_position",
            "time_source": exp.time_source, "time_assumptions": exp.time_assumptions}


def spectrum(exps: list[Exposure], coord_at: Callable[[float], SkyCoord]) -> list[dict]:
    """Точки спектра со всех экспозиций. coord_at(mjd) учитывает собственное движение цели."""
    pts = [p for e in exps if (p := forced_photometry(e, coord_at(e.mjd))) is not None]
    return sorted(pts, key=lambda p: p["wave_um"] if p["wave_um"] is not None else float("inf"))


def light_curve(points: list[dict], wave_min: float, wave_max: float) -> list[dict]:
    """Поток по времени в узком окне λ: иначе переменность спутается с формой спектра (LVF)."""
    sel = [p for p in points if p["wave_um"] is not None and wave_min <= p["wave_um"] <= wave_max and p["flux_ujy"] is not None]
    return sorted(sel, key=lambda p: p["mjd"])
