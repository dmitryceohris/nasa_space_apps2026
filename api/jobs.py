"""Background jobs, one per cache key. The front polls the same URL until the job is ready (202 -> 200).

Image processing is NOT done here: the frames come from src/processing (coadd with zodiacal-light subtraction,
bad-pixel masks, wavelength at the target, difference, trail) and the spectrum from build_spectrum. This module
only finds the exposures (local table or IRSA TAP), downloads the cutouts, calls those functions, writes the
PNG files and adds the object centres.
"""
import io
import hashlib
import json
import logging
import os
import threading
import time
import uuid
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from Exposure import DETECTOR_BANDS, load_exposure
from processing import render as prender
from processing.epochs import epoch_limits, select_per_epoch, smear_px
from processing.frames import (BACKGROUND_METHOD, coadd, difference_details, grid_header_json, make_grid,
                               spectral_sampling, trail)
from processing.pipeline import CREDITS, build_spectrum, lightcurve_from_spectrum

from . import cache, centroid, irsa, render, resolve
from .errors import ApiError, PUBLIC_INTERNAL_MESSAGE
from .params import Params

log = logging.getLogger("jobs")
POOL = ThreadPoolExecutor(max_workers=4)
JOBS = {}
LOCK = threading.Lock()
SMEAR_WARN_PX = 0.5  # known proper motion moves the object more than this within one epoch -> warning
ASSOC_PX = 1.0  # Conservative display gate for candidate disagreement; NOT a calibrated association radius.
# A match within this radius does not establish identity or registration accuracy, especially for blends
# and fields without usable reference stars. All candidates remain available separately.
warnings.filterwarnings("ignore")


@dataclass
class Job:
    key: str
    id: str
    status: str = "running"  # running | error (every finished result goes to the cache, empty ones with a short TTL)
    stage: str = "query"  # query -> download -> process -> render
    done: int = 0
    total: int = 0
    message: str = ""
    created: float = field(default_factory=time.monotonic)
    started: bool = False
    finished: float | None = None


def _poll(k: str, runner):
    """-> (http_status, payload): 200 cached/ready/error, 202 running. Starts runner(job) if nothing runs for k."""
    res = cache.read(k)
    if res:
        return 200, res
    with LOCK:
        # The worker may have published and removed its Job since the first cache read.
        res = cache.read(k)
        if res:
            return 200, res
        now = time.monotonic()
        retention = float(os.environ.get("JOB_RETENTION_SECONDS", "300"))
        for old_key, old in list(JOBS.items()):
            if old.status != "running" and now - (old.finished or old.created) > retention:
                JOBS.pop(old_key)
        j = JOBS.get(k)
        if j is None:
            # Only queued/running work counts: retained failures wait to be reported, they hold no worker.
            if sum(o.status == "running" for o in JOBS.values()) >= int(os.environ.get("MAX_PENDING_JOBS", "16")):
                raise ApiError(503, "busy", "Очередь запросов заполнена. Повторите запрос позже")
            j = JOBS[k] = Job(k, uuid.uuid4().hex[:8])
            POOL.submit(_run, j, runner)
        if j.status == "error":  # drop it so the next request retries
            JOBS.pop(k)
            return 200, {"status": "error", "message": j.message}
        if now - j.created > float(os.environ.get("MAX_JOB_SECONDS", "1800")):
            raise ApiError(503, "timeout", "Задача превысила время ожидания; загрузки ещё завершаются. Повторите запрос позже")
        return 202, {"status": "running", "job_id": j.id, "stage": j.stage if j.started else "queued", "done": j.done, "total": j.total,
                     "message": j.message}


def _run(j, runner):
    try:
        j.started = True
        _deadline(j)
        res = runner(j)
        _deadline(j)
        cache.write(j.key, res)
        with LOCK:
            if JOBS.get(j.key) is j:
                JOBS.pop(j.key)
    except ApiError as e:
        with LOCK:
            j.status, j.message, j.finished = "error", e.message, time.monotonic()
    except Exception as e:
        log.exception("job %s failed", j.key)
        with LOCK:
            j.status, j.message, j.finished = "error", PUBLIC_INTERNAL_MESSAGE, time.monotonic()


def _stage(j, stage, msg="", done=0, total=0):
    _deadline(j)
    j.stage, j.message, j.done, j.total = stage, msg, done, total


def _deadline(j):
    if time.monotonic() - j.created > float(os.environ.get("MAX_JOB_SECONDS", "1800")):
        raise ApiError(503, "timeout", "Превышено время обработки запроса")


# ---------------------------------------------------------------- view
def _model_key(key, r):
    # Fixed fields and tracking queries at the same centre must never share results.
    if not r.get("motion"):
        return key
    model = json.dumps(r["motion"], sort_keys=True, ensure_ascii=False)
    return hashlib.sha1((key + model).encode()).hexdigest()[:16]


def _view_key(p: Params, r: dict) -> str:
    ra, dec = r["center"]
    key = cache.key(ra, dec, p.fov, p.det, p.mode, p.max_per_epoch)
    key = hashlib.sha1((key + irsa.selection_token()).encode()).hexdigest()[:16]
    return _model_key(key, r)


def view(p: Params):
    r = resolve.for_params(p)
    return _poll(_view_key(p, r), lambda j: _pipeline(j, p, r))


def cached_result(p: Params):
    return cache.read(_view_key(p, resolve.for_params(p)))


def gif(p: Params) -> bytes:
    res = cached_result(p)
    if not res or not res.get("frames"):
        raise ApiError(409, "not_ready", "результат ещё не готов: сначала откройте /api/view")
    paths = [cache.CACHE / f["image_url"].removeprefix("/cache/") for f in res["frames"]]
    buf = io.BytesIO()
    render.save_gif(paths, buf)
    return buf.getvalue()


def prewarm():
    """Start the 5 presets with their own det/fov (mode coadd) so a demo never downloads live."""
    for t in resolve.PRESETS:
        try:
            view(Params(t["id"], None, None, float(t["fov"]), t["det"], "coadd", 6))
        except Exception:
            log.exception("prewarm %s failed", t["id"])


def _motion(r):
    """Reference proper motion from SIMBAD. Measured motion is not computed yet: it stays null."""
    o = r["motion"]
    if not o or o.get("pmra_masyr") is None or o.get("pmdec_masyr") is None:
        return None
    return {"pm_reference_arcsec_yr": [round(o["pmra_masyr"] / 1000, 4), round(o["pmdec_masyr"] / 1000, 4)],
            "source": o["source"], "pm_measured_arcsec_yr": None}


def _pm_total(r):
    """Total reference proper motion, arcsec/yr, or None."""
    o = r["motion"]
    if not o or o.get("pmra_masyr") is None or o.get("pmdec_masyr") is None:
        return None
    return float(np.hypot(o["pmra_masyr"], o["pmdec_masyr"])) / 1000


def _processing(p):
    """How the frames and centres were obtained, shown to the user next to the result."""
    from processing.frames import DIFF_METHOD

    return {"epochs": epoch_limits(), "background": BACKGROUND_METHOD,
            "registration": {"name": "header_wcs_reprojection", "independently_verified": False,
                             "astrometric_precision_arcsec": None},
            "difference": DIFF_METHOD if p.mode in ("diff", "trail") else None,
            "measurement": centroid.METHOD}


def _sky(r, mjd):
    """Position of the object on `mjd`: proper motion if known, else the fixed query centre."""
    if r["motion"]:
        ra, dec = resolve.propagate(r["motion"], mjd)
        return SkyCoord(ra * u.deg, dec * u.deg)
    return SkyCoord(r["center"][0] * u.deg, r["center"][1] * u.deg)


def _pipeline(j, p, r):
    ra, dec = r["center"]
    center = SkyCoord(ra * u.deg, dec * u.deg)
    wcs, shape = make_grid(center, p.fov)
    res = {"status": "ready",
           "query": {"ra": ra, "dec": dec, "fov": p.fov, "det": p.det, "mode": p.mode, "target": r["label"]},
           "grid": grid_header_json(wcs, shape), "frames": [], "trail_url": None, "motion": _motion(r),
           "target_info": r.get("motion"), "center_mjd": r.get("center_mjd"),
           "data_status": "no_archive_results", "warnings": [], "credits": CREDITS, "processing": _processing(p)}

    _stage(j, "query", "Поиск снимков")
    meta = irsa.exposures(ra, dec, p.det)
    res["archive_source"] = getattr(meta, "source", "tables" if irsa.selection_token() != "archive" else "tap")
    res["warnings"].extend(getattr(meta, "warnings", []))
    if not meta:
        res["warnings"].append("По этим параметрам нет доступных кадров в ответе архива")
        return res
    groups = select_per_epoch([m["mjd"] for m in meta], p.max_per_epoch)  # epochs, <= max_per_epoch each
    if p.mode == "single":  # one exposure per epoch: download only those
        groups = [g[:1] for g in groups]
    flat = [i for g in groups for i in g]

    size = p.fov / 60 * 1.5  # cutout a bit larger than the grid so rotated footprints still fill it
    _stage(j, "download", f"Скачано 0 из {len(flat)} снимков", 0, len(flat))
    paths, failed = {}, 0
    ex = ThreadPoolExecutor(6)
    try:
        futs = {ex.submit(irsa.cutout, meta[i]["uri"], ra, dec, size): i for i in flat}
        for f in as_completed(futs):
            _deadline(j)
            try:
                paths[futs[f]] = f.result()
            except Exception as e:
                failed += 1
                log.warning("cutout failed: %s", e)
            j.done += 1
            j.message = f"Скачано {j.done} из {j.total} снимков"
    finally:
        # Do not drain the remaining queue after the request deadline. In-flight HTTP calls are bounded.
        ex.shutdown(wait=True, cancel_futures=True)

    _stage(j, "process", "Обработка", 0, len(groups))
    frames = []
    for g in groups:
        _deadline(j)
        exps = []
        for i in g:
            if i not in paths:
                continue
            try:
                exps.append(load_exposure(str(paths[i]), mjd_hint=meta[i]["mjd"],
                                          obs_id_hint=irsa.obs_id(meta[i]["uri"]) or ""))
            except Exception as e:
                failed += 1
                log.warning("load failed: %s", e)
        if exps:
            # lambda is read at the target position on the epoch date, not at the field centre
            f = coadd(exps, wcs, shape, _sky(r, float(np.mean([e.mjd for e in exps]))))
            if f is not None:
                frames.append(f)
        j.done += 1
    if failed:
        res["warnings"].append(f"Не удалось получить или обработать {failed} из {len(flat)} снимков")

    res["data_status"] = ("download_failed" if not paths else "unusable_data") if not frames else (
        "partial" if failed or sum(f.n_exp for f in frames) < len(flat) else "single_epoch" if len(frames) == 1 else "available")
    used = sum(f.n_exp for f in frames)
    res["exposure_counts"] = {"selected": len(flat), "failed": failed, "used": used,
                              "rejected": len(flat) - failed - used}
    if not frames:
        res["warnings"].append("Все загрузки выбранных экспозиций не удались" if not paths else
                               "Файлы получены, но пригодных для изображения пикселей нет")
    if res["exposure_counts"]["rejected"]:
        res["warnings"].append(f"Отброшено непригодных экспозиций: {res['exposure_counts']['rejected']}")
    if len(frames) == 1:
        res["warnings"].append("Доступна одна эпоха: последовательность разных дат пока недоступна")
    if frames:
        res["warnings"].extend([
            "Независимая регистрация по опорным звёздам не проверена; точность измеренного движения не подтверждена",
            "Центроиды — кандидаты, не подтверждённые положения цели; ошибки центроидов не включают все систематические эффекты",
        ])
    if len(frames) > 1:
        # Even a <5% change can select a different narrow LVF band. The average wavelength
        # and detector ID cannot certify matching bandpasses, PSFs or photometric response.
        res["warnings"].append(
            "Различия яркости и разностные остатки не подтверждают переменность или движение: "
            "нужны сопоставимые полосы/отклики, проверка масок и точности совмещения")
    diffs = template = scales = None
    if p.mode in ("diff", "trail"):
        if len(frames) < 2:
            res["warnings"].append(f"{p.mode} требует минимум 2 эпохи")
            if p.mode == "diff":  # nothing meaningful to show
                frames = []
        else:
            dd = difference_details(frames)
            diffs, template, scales = dd["diffs"], dd["template"], dd["noise_scale"]

    _stage(j, "render", "Рисую кадры", 0, len(frames))
    d = cache.CACHE / j.key
    d.mkdir(exist_ok=True)
    vmin, vmax = prender.shared_limits([f.image for f in frames]) if frames else (0.0, 1.0)
    h = shape[0]
    for i, f in enumerate(frames):
        signed = p.mode == "diff"
        arr = diffs[i] if signed else f.image
        if not np.isfinite(arr).any():
            res["warnings"].append(f"Эпоха {i + 1}: нет пригодных пикселей для " + ("разности" if signed else "изображения"))
            continue
        if signed:
            prender.save_diff(arr, d / f"f{i}.png")
        else:
            prender.save_gray(arr, d / f"f{i}.png", vmin, vmax, saturated=f.contam)
        prender.save_mask(f.contam, d / f"f{i}_flags.png")  # bad and contaminated pixels; never SOURCE
        pred = pred_fits = None
        if r["motion"]:  # where the object should be on the frame date
            x, y = wcs.world_to_pixel(_sky(r, f.mjd))
            if np.isfinite([x, y]).all() and 0 <= x < shape[1] and 0 <= y < h:
                pred = prender.fits_to_png_xy(x, y, h)
                pred_fits = (float(x), float(y))
            else:
                res["warnings"].append(f"Эпоха {i + 1}: прогноз вне поля; маркер скрыт")
        m = centroid.measure(arr, signed, saturated=f.contam)
        # Suppress a clearly discrepant candidate in the legacy single-marker field. Agreement alone does
        # not identify the target (WISE 0855 at 4.42 um had a candidate 17 px from its prediction).
        if not r["motion"]:
            assoc = "no_prediction"
        elif pred_fits is None:  # the model puts the target outside this frame: no object here is the target
            assoc = "prediction_outside_field"
        else:
            assoc = "no_measurement"
        offset = None
        if m["xy"] is not None and pred_fits is not None:
            offset = float(np.hypot(m["xy"][0] - pred_fits[0], m["xy"][1] - pred_fits[1]))
            assoc = "far_from_prediction" if offset > ASSOC_PX else "consistent_with_prediction"
        cands = [{"kind": c["kind"], "sign": c["sign"], "xy": prender.fits_to_png_xy(c["xy"][0], c["xy"][1], h),
                  "err_px": [round(e, 3) for e in c["err_px"]], "area_px": c["area"], "snr": round(c["snr"], 1),
                  "quality": c["quality"], "flags": c["flags"]}
                 for c in m["candidates"]]
        coverage = float(np.isfinite(arr).sum() / arr.size)
        limitations = []
        if coverage < 1:
            limitations.append("incomplete_pixel_coverage")
        if f.contam.any():
            limitations.append("contaminated_pixels")
        if any(s["source"] != "spectral_wcs" for s in f.wavelength_samples):
            limitations.append("approximate_wavelength")
        if any(s["wave_um"] is None or s["bw_um"] is None for s in f.wavelength_samples):
            limitations.append("invalid_wavelength")
        if any(t["assumptions"] for t in f.time_sources):
            limitations.append("time_assumptions")
        res["frames"].append({
            "index": i, "label": f"Эпоха {i + 1}", "mjd": round(f.mjd, 4), "date": f.date,
            "mjd_span": [round(f.mjd_min, 4), round(f.mjd_max, 4)], "det": f.detector,
            "wave_um": round(f.wave_um, 4) if np.isfinite(f.wave_um) else None, "n_exp": f.n_exp,
            "image_url": f"/cache/{j.key}/f{i}.png", "flags_url": f"/cache/{j.key}/f{i}_flags.png",
            "predicted_px": pred,
            "prediction_status": ("in_frame" if pred else "outside_field") if r["motion"] else "no_model",
            # measured_px: the centre when the choice is unambiguous (status ok); null when nothing is found or
            # when the largest and the central object differ (both are listed, nothing is guessed). In diff mode
            # only the "+" residual (where the object is in this frame) gives the position.
            "measured_px": prender.fits_to_png_xy(m["xy"][0], m["xy"][1], h)
            if m["status"] == "ok" and m["quality"] != "poor"
            and assoc not in ("far_from_prediction", "prediction_outside_field") else None,
            "measured_candidates": cands,
            "measurement": {"status": m["status"], "reason": m["reason"], "quality": m["quality"],
                            "flags": m["flags"], "n_rejected_hot": m["n_rejected_hot"],
                            "n_low_snr": m["n_low_snr"], "target_association": assoc,
                            "offset_from_prediction_px": None if offset is None else round(offset, 3),
                            "err_px": None if m["err_px"] is None else [round(e, 3) for e in m["err_px"]]},
            "quality": {"status": "limited" if limitations else "unchecked", "verified": False,
                        "coverage_fraction": round(coverage, 5), "reasons": limitations},
            "time_sources": f.time_sources, "wavelength_samples": f.wavelength_samples,
            "span_days": round(f.span_days, 3), "bg_offsets": f.bg_offsets,
            "diff_noise_scale": scales[i] if scales else None,
            "obs_ids": f.obs_ids})
    if frames and not res["frames"]:
        res["data_status"] = "unusable_data"
    elif len(res["frames"]) < len(frames):
        res["data_status"] = "partial"
    if diffs is not None and res["frames"]:
        prender.save_rgb(trail(diffs), d / "trail.png", background=template, vmin=vmin, vmax=vmax)
        res["trail_url"] = f"/cache/{j.key}/trail.png"

    ss = res["spectral_sampling"] = spectral_sampling(frames)
    if ss["different_channels"]:
        width = f", ширина канала ≈ {ss['channel_width_um']:.3f} мкм" if ss["channel_width_um"] else ""
        res["warnings"].append(f"Цель снята в разных спектральных каналах ({ss['wave_min_um']:.2f}-{ss['wave_max_um']:.2f} мкм"
                               f"{width}): разница яркости может быть связана с формой спектра; "
                               "переменность и движение по этим остаткам не подтверждены")
    c = (h // 2, shape[1] // 2)
    if any(f.contam[c[0] - 1:c[0] + 2, c[1] - 1:c[1] + 2].any() for f in frames):
        res["warnings"].append("В центре поля есть насыщенные или загрязнённые пиксели (слой флагов)")
    pm = _pm_total(r)
    if pm and frames:
        pix = res["grid"]["pixel_scale_arcsec"]
        worst = max(frames, key=lambda f: f.span_days)
        sm = smear_px(pm, worst.span_days, pix)
        if sm > SMEAR_WARN_PX:
            res["warnings"].append(f"Объект смещается на {sm:.1f} пикс. за эпоху ({worst.span_days:.0f} сут): "
                                   "внутри кадра он размазан, центр — среднее положение за эпоху")
    return res


# ---------------------------------------------------------------- spectrum and light curve
def _spectrum_job(j, r):
    ra, dec = r["center"]
    center = SkyCoord(ra * u.deg, dec * u.deg)

    def progress(stage, done, total, msg):
        _stage(j, stage, msg, done, total)

    from Exposure.irsa import ArchiveQueryError

    try:
        res = build_spectrum(center, lambda mjd: _sky(r, mjd), progress=progress)
    except ArchiveQueryError as error:
        raise ApiError(502, "upstream", str(error)) from error
    res["query"] = {"ra": ra, "dec": dec, "target": r["label"]}
    res["target_info"] = r.get("motion")
    res["center_mjd"] = r.get("center_mjd")
    return res


def spectrum_view(p: Params):
    """Only target / ra+dec of p are used."""
    r = resolve.for_params(p)
    return _poll(_model_key(cache.key_spectrum(*r["center"]), r), lambda j: _spectrum_job(j, r))


def default_window(det: int):
    """+-5% around the middle of the detector band (spec 3.6)."""
    lo, hi, _ = DETECTOR_BANDS[det]
    mid = (lo + hi) / 2
    return mid * 0.95, mid * 1.05


def lightcurve_view(p: Params, wave_min: float, wave_max: float):
    """Light curve = a window of the cached spectrum; while the spectrum job runs the client gets its 202."""
    status, payload = spectrum_view(p)
    if status != 200 or payload.get("status") != "ready":
        return status, payload
    res = lightcurve_from_spectrum(payload, wave_min, wave_max)
    res["query"] = payload.get("query")
    res["target_info"] = payload.get("target_info")
    res["warnings"] = list(payload.get("warnings", []))  # copy: never mutate the cached spectrum
    res["warnings"] += res.pop("lightcurve_warnings", [])  # points from different spectral channels
    if not res["flux"]:
        res["warnings"].append(f"В окне {wave_min:.3f}–{wave_max:.3f} мкм нет надёжных точек спектра; "
                               "задайте wave_min и wave_max шире или на другом детекторе")
    return 200, res
