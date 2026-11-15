"""Полная цепочка: координаты -> снимки IRSA -> эпохи -> кадры -> PNG + JSON для фронтенда."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord

from Exposure import DETECTOR_BANDS, PIXEL_SCALE_ARCSEC, Exposure, download_rows, query_spherex

from . import render
from .epochs import epoch_limits, group_epochs, select_per_epoch
from .frames import (BACKGROUND_METHOD, DIFF_METHOD, Frame, coadd, difference, grid_header_json, make_grid,
                     spectral_sampling, trail)
from .photometry import PHOTOMETRY_METHOD, light_curve, spectrum
from .targets import Target
from .validation import accepted

Progress = Callable[[str, int, int, str], None]  # (stage, done, total, message)
CREDITS = "NASA SPHEREx via NASA/IPAC IRSA"


def _noop(*_a, **_k):
    pass


# ------------------------------------------------------------------ кадры («видео»)
def frames_from_exposures(exps: list[Exposure], center: SkyCoord, fov_arcmin: float,
                          gap_days: float = 30, max_span_days: float = 20):
    wcs, shape = make_grid(center, fov_arcmin)
    groups = group_epochs([e.mjd for e in exps], gap_days, max_span_days)
    frames = [f for g in groups if (f := coadd([exps[i] for i in g], wcs, shape, center)) is not None]
    return frames, wcs, shape


def render_view(frames: list[Frame], wcs, shape, out_dir: Path, url_prefix: str,
                target: Target | None = None) -> dict:
    """Пишет PNG всех режимов и возвращает JSON-ответ /api/view (раздел 3.4 docs/API.md)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    h = shape[0]
    warnings: list[str] = []
    if not frames:
        return {"status": "ready", "grid": grid_header_json(wcs, shape), "frames": [],
                "trail_url": None, "gif_url": None, "motion": None, "credits": CREDITS,
                "warnings": ["SPHEREx ещё не наблюдал эту область на этом детекторе"]}

    vmin, vmax = render.shared_limits([f.image for f in frames])
    diffs, template = difference(frames) if len(frames) > 1 else ([np.zeros(shape)], frames[0].image)

    out_frames, coadd_pngs = [], []
    for i, (f, d) in enumerate(zip(frames, diffs)):
        p_img, p_diff, p_flags = out_dir / f"f{i}.png", out_dir / f"f{i}_diff.png", out_dir / f"f{i}_flags.png"
        render.save_gray(f.image, p_img, vmin, vmax, saturated=f.contam)
        render.save_diff(d, p_diff)
        render.save_mask(f.contam, p_flags)
        coadd_pngs.append(p_img)
        pred = None
        if target is not None:
            x, y = wcs.world_to_pixel(target.coord_at(f.mjd))
            pred = render.fits_to_png_xy(x, y, h)
        out_frames.append({
            "index": i, "label": f"Эпоха {i + 1}", "mjd": round(f.mjd, 4), "date": f.date,
            "mjd_span": [round(f.mjd_min, 4), round(f.mjd_max, 4)],
            "det": f.detector, "wave_um": None if not np.isfinite(f.wave_um) else round(f.wave_um, 4),
            "n_exp": f.n_exp, "obs_ids": f.obs_ids,
            "urls": {"coadd": f"{url_prefix}/{p_img.name}", "diff": f"{url_prefix}/{p_diff.name}",
                     "flags": f"{url_prefix}/{p_flags.name}"},
            "predicted_px": pred, "measured_px": None,
        })

    trail_png = out_dir / "trail.png"
    render.save_rgb(trail(diffs), trail_png, background=template, vmin=vmin, vmax=vmax)
    gif = out_dir / "movie.gif"
    render.save_gif(coadd_pngs, gif)

    # предупреждения, которые честно объясняют пользователю ограничения данных
    ss = spectral_sampling(frames)
    if ss["different_channels"]:
        warnings.append(f"Цель снята в разных спектральных каналах ({ss['wave_min_um']:.2f}–{ss['wave_max_um']:.2f} мкм): "
                        "разница яркости может быть связана с формой спектра; переменность и движение не подтверждены")
    if len(frames) < 2:
        warnings.append("Только одна эпоха — движение пока не видно")
    c = (shape[0] // 2, shape[1] // 2)
    if any(f.contam[c[0] - 1:c[0] + 2, c[1] - 1:c[1] + 2].any() for f in frames):
        warnings.append("В центре поля есть насыщенные или загрязнённые пиксели (слой «флаги»)")

    motion = None
    if target is not None and (target.pmra or target.pmdec):
        motion = {"pm_reference_arcsec_yr": [round(target.pmra / 1000, 4), round(target.pmdec / 1000, 4)],
                  "source": target.source, "pm_measured_arcsec_yr": None}
    return {"status": "ready", "grid": grid_header_json(wcs, shape), "frames": out_frames,
            "trail_url": f"{url_prefix}/{trail_png.name}", "gif_url": f"{url_prefix}/{gif.name}",
            "motion": motion, "warnings": warnings, "credits": CREDITS,
            "processing": {"epochs": epoch_limits(), "background": BACKGROUND_METHOD,
                           "difference": DIFF_METHOD if len(frames) > 1 else None}}


def build_view(center: SkyCoord, fov_arcmin: float, det: int, max_per_epoch: int,
               out_dir: Path, url_prefix: str, target: Target | None = None,
               progress: Progress = _noop) -> dict:
    progress("query", 0, 0, "Ищем снимки SPHEREx в IRSA…")
    tbl = query_spherex(center, detector=det)
    if len(tbl) == 0:
        wcs, shape = make_grid(center, fov_arcmin)
        return render_view([], wcs, shape, out_dir, url_prefix, target)
    sel = [i for g in select_per_epoch(np.asarray(tbl["mjd"], float), max_per_epoch) for i in g]
    rows = [tbl[i] for i in sel]
    size_deg = fov_arcmin / 60 * 1.45 + 0.01  # с запасом: квадрат сетки повёрнут относительно кадра
    exps = download_rows(rows, center, size_deg,
                         progress=lambda d, t: progress("download", d, t, f"Скачано {d} из {t} снимков"))
    progress("process", 0, 0, "Складываем эпохи и строим кадры…")
    frames, wcs, shape = frames_from_exposures(exps, center, fov_arcmin)
    progress("render", 0, 0, "Рисуем картинки…")
    res = render_view(frames, wcs, shape, out_dir, url_prefix, target)
    res["available_exposures"] = int(len(tbl))
    return res


# ------------------------------------------------------------------ спектр и кривая блеска
def build_spectrum(center: SkyCoord, coord_at: Callable[[float], SkyCoord], per_epoch: int = 2,
                   progress: Progress = _noop) -> dict:
    """Спектр по всем 6 детекторам: маленькие вырезки, per_epoch снимков на эпоху и детектор."""
    progress("query", 0, 0, "Ищем снимки всех детекторов…")
    tbl = query_spherex(center)
    rows = []
    for det in sorted(set(int(d) for d in tbl["detector"])):
        sub = tbl[tbl["detector"] == det]
        for g in select_per_epoch(np.asarray(sub["mjd"], float), per_epoch):
            rows += [sub[i] for i in g]
    max_offset = max((float(center.separation(coord_at(float(row["mjd"]))).arcsec) for row in rows), default=0.0)
    # Include the complete background annulus at the propagated location on every selected date.
    size_deg = max(0.035, 2 * (max_offset + PHOTOMETRY_METHOD["params"]["r_out_px"] * PIXEL_SCALE_ARCSEC) * 1.1 / 3600)
    exps = download_rows(rows, center, size_deg,
                         progress=lambda d, t: progress("download", d, t, f"Скачано {d} из {t} снимков"))
    pts = spectrum(exps, coord_at)
    good = _good(pts)
    warnings = []
    if len(exps) < len(rows):
        warnings.append(f"Не удалось загрузить {len(rows) - len(exps)} из {len(rows)} экспозиций")
    if len(pts) < len(exps):
        warnings.append(f"Для {len(exps) - len(pts)} экспозиций нет пригодного покрытия апертуры/фонового кольца")
    if any(p["wavelength_source"] not in ("spectral_wcs", "cwave_cband") for p in pts):
        warnings.append("Длина волны части точек приблизительна: спектральный WCS недоступен")
    if not rows:
        warnings.append("По этим параметрам нет доступных кадров в ответе архива")
    elif not exps:
        warnings.append("Все загрузки выбранных экспозиций не удались")
    elif not good:
        warnings.append("Надёжных измерений для графика нет; проверьте покрытие, маски и флаги качества")
    state = ("no_archive_results" if not rows else "download_failed" if not exps else
             "unusable_data" if not good else "partial" if len(exps) < len(rows) else "available")
    counts = {s:sum(p["validation"]["status"]==s for p in pts) for s in ("validated","provisional","rejected")}
    return {"status": "ready", "points": pts, "n_points": len(pts), "warnings": warnings,
            "quality_counts": counts,
            "data_status": state, "exposure_counts": {"selected": len(rows), "downloaded": len(exps),
                                                        "with_measurements": len(pts), "reliable": len(good)},
            "cutout_size_deg": size_deg, "max_prediction_offset_arcsec": max_offset,
            "n_saturated": sum(1 for p in pts if {"OVERFLOW", "NONLINEAR"} & set(p["flags"])),
            "n_flagged": sum(1 for p in pts if not p["ok"]), "credits": CREDITS,
            "method": {"name":"native_epsf_with_preliminary_aperture_fallback",
                       "description":"Native ePSF measurements when calibrated; aperture diagnostics remain provisional",
                       "limitations":PHOTOMETRY_METHOD["limitations"],"not_implemented":[]},
            "bands": {str(k): [v[0], v[1]] for k, v in DETECTOR_BANDS.items()},
            # плоские массивы только надёжных точек (ok): их напрямую рисует frontend/index.html
            "wavelength": [p["wave_um"] for p in good], "flux": [p["flux_ujy"] for p in good],
            "flux_err": [p["err_ujy"] for p in good]}


def _good(pts: list[dict]) -> list[dict]:
    return [p for p in pts if accepted(p)]


def lightcurve_from_spectrum(spec: dict, wave_min: float, wave_max: float) -> dict:
    from astropy.time import Time

    diagnostics = light_curve(spec["points"], wave_min, wave_max, include_diagnostics=True)
    pts = light_curve(spec["points"], wave_min, wave_max)
    good = _good(pts)
    w = np.array([p["wave_um"] for p in diagnostics if p.get("wave_um") is not None], float)
    bw = np.array([p["bw_um"] for p in diagnostics if p.get("bw_um")], float)
    width = float(np.median(bw)) if bw.size else None
    mixed = bool(w.size > 1 and width and np.ptp(w) > width)
    lc_warn = []
    if mixed:
        lc_warn.append(f"Точки кривой блеска сняты в разных спектральных каналах ({w.min():.3f}–{w.max():.3f} мкм при "
                       f"ширине канала ≈ {width:.3f} мкм): изменение потока не доказывает переменность")
    from .spectral import compatible_groups
    groups = compatible_groups([p for p in diagnostics if accepted(p)])
    series = [{"dates":[Time(p["mjd"],format="mjd").isot for p in g],
               "flux":[p["flux_ujy"] for p in g],"flux_err":[p["err_ujy"] for p in g],
               "obs_ids":[p["obs_id"] for p in g]} for g in groups]
    if diagnostics and not groups:
        lc_warn.append("Нет принятых измерений с подтверждённой сопоставимостью спектральных откликов")
    return {"status": "ready", "wave_window_um": [wave_min, wave_max], "points": pts,
            "diagnostic_points": diagnostics,
            "series":series,"spectral_comparison":{"status":"validated" if groups else "provisional",
                                                       "n_groups":len(groups)},
            "spectral_sampling": {"wave_min_um": float(w.min()) if w.size else None,
                                  "wave_max_um": float(w.max()) if w.size else None,
                                  "channel_width_um": width, "different_channels": mixed},
            "lightcurve_warnings": lc_warn,
            "method": dict(PHOTOMETRY_METHOD, window="поток по времени только в окне λ: иначе переменность "
                                                     "спуталась бы с формой спектра (LVF)"),
            "dates": [Time(p["mjd"], format="mjd").to_value("iso", subfmt="date") for p in good],
            "flux": [p["flux_ujy"] for p in good], "flux_err": [p["err_ujy"] for p in good]}
