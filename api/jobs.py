"""Background jobs, one per cache key. The front polls the same URL until the job is ready (202 -> 200).

Image processing is NOT done here: the frames come from src/processing (coadd with zodiacal-light subtraction,
bad-pixel masks, wavelength at the target, difference, trail) and the spectrum from build_spectrum. This module
only finds the exposures (local table or IRSA TAP), downloads the cutouts, calls those functions, writes the
PNG files and adds the object centres.
"""
import io
import logging
import os
import threading
import time
import uuid
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from Exposure import DETECTOR_BANDS, load_exposure
from processing import render as prender
from processing.epochs import epoch_limits, select_per_epoch, smear_px
from processing.frames import BACKGROUND_METHOD, coadd, difference_details, grid_header_json, make_grid, trail
from processing.pipeline import CREDITS, build_spectrum, lightcurve_from_spectrum

from . import cache, centroid, irsa, render, resolve
from .errors import ApiError
from .params import Params

log = logging.getLogger("jobs")
POOL = ThreadPoolExecutor(max_workers=4)
JOBS = {}
LOCK = threading.Lock()
_active_jobs = 0
WAVE_SPREAD_WARN = 0.05  # same threshold as src/processing/pipeline.py: lambda change between epochs
SMEAR_WARN_PX = 0.5  # known proper motion moves the object more than this within one epoch -> warning
warnings.filterwarnings("ignore")


@dataclass
class Job:
    key: str
    id: str
    status: str = "running"  # running | ready (transient/partial/stale result) | error
    stage: str = "query"  # query -> download -> process -> render
    done: int = 0
    total: int = 0
    message: str = ""
    result: dict | None = None
    hold_until: float = 0.0


def _max_active() -> int:
    return int(os.environ.get("JOB_MAX_ACTIVE", "16"))


def _hold_sec() -> float:
    return float(os.environ.get("JOB_HOLD_SEC", "2.0"))


def _empty_ttl() -> float:
    return float(os.environ.get("CACHE_EMPTY_TTL", "300.0"))


def _retry_after() -> str:
    return os.environ.get("JOB_RETRY_AFTER", "5")


def _poll(k: str, runner):
    """-> (http_status, payload): 200 cached/ready/error, 202 running. Starts runner(job) if nothing runs for k."""
    global _active_jobs
    with LOCK:
        res = cache.read(k)
        if res is not None:
            JOBS.pop(k, None)
            return 200, res
        now = time.monotonic()
        j = JOBS.get(k)
        if j is not None and j.status in ("error", "ready") and now >= j.hold_until:
            JOBS.pop(k, None)
            j = None
        if j is None:
            if _active_jobs >= _max_active():
                raise ApiError(429, "queue_full", "очередь задач переполнена, повторите запрос позже",
                               headers={"Retry-After": _retry_after()})
            _active_jobs += 1
            j = JOBS[k] = Job(k, uuid.uuid4().hex[:8])
            try:
                POOL.submit(_run, j, runner)
            except Exception:
                _active_jobs -= 1
                JOBS.pop(k, None)
                raise
        if j.status == "error":
            return 200, {"status": "error", "message": j.message}
        if j.status == "ready":
            return 200, j.result
        return 202, {"status": "running", "job_id": j.id, "stage": j.stage, "done": j.done, "total": j.total,
                     "message": j.message}


def _stale_fallback(k: str, msg: str) -> dict | None:
    stale = cache.read(k, allow_stale=True)
    if stale is None:
        return None
    stale = dict(stale)
    stale["stale"] = True
    warns = list(stale.get("warnings") or [])
    w_msg = f"Архив недоступен ({msg}); показан устаревший результат из кэша"
    if w_msg not in warns:
        warns.append(w_msg)
    stale["warnings"] = warns
    return stale


def _run(j, runner):
    global _active_jobs
    try:
        res = runner(j)
        partial = bool(res.pop("_partial", False))
        empty = (res.get("frames") is not None and not res["frames"]) or (
            res.get("points") is not None and not res["points"])
        if partial:
            with LOCK:
                j.result, j.status, j.hold_until = res, "ready", time.monotonic() + _hold_sec()
        elif empty:
            cache.write(j.key, res, ttl=_empty_ttl())
            with LOCK:
                JOBS.pop(j.key, None)
        else:
            cache.write(j.key, res)
            with LOCK:
                JOBS.pop(j.key, None)
    except ApiError as e:
        if e.status == 502:
            stale = _stale_fallback(j.key, e.message)
            if stale is not None:
                with LOCK:
                    j.result, j.status, j.hold_until = stale, "ready", time.monotonic() + _hold_sec()
                return
        with LOCK:
            j.status, j.message, j.hold_until = "error", e.message, time.monotonic() + _hold_sec()
    except Exception as e:
        log.exception("job %s failed", j.key)
        with LOCK:
            j.status, j.message, j.hold_until = "error", f"внутренняя ошибка: {e}", time.monotonic() + _hold_sec()
    finally:
        with LOCK:
            _active_jobs = max(0, _active_jobs - 1)


def _stage(j, stage, msg="", done=0, total=0):
    j.stage, j.message, j.done, j.total = stage, msg, done, total


def _atomic_png(save_fn, target_path, first_arg, *args, **kwargs):
    target_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = target_path.parent / f".tmp_{uuid.uuid4().hex[:8]}_{target_path.name}"
    save_fn(first_arg, tmp, *args, **kwargs)
    tmp.replace(target_path)


# ---------------------------------------------------------------- view
def _view_key(p: Params, r: dict) -> str:
    ra, dec = r["center"]
    return cache.key(ra, dec, p.fov, p.det, p.mode, p.max_per_epoch)


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
           "warnings": [], "credits": CREDITS, "processing": _processing(p)}

    _stage(j, "query", "Поиск снимков")
    meta = irsa.exposures(ra, dec, p.det)
    if not meta:
        res["warnings"].append("SPHEREx ещё не наблюдал эту область")
        return res
    groups = select_per_epoch([m["mjd"] for m in meta], p.max_per_epoch)  # epochs, <= max_per_epoch each
    if p.mode == "single":  # one exposure per epoch: download only those
        groups = [g[:1] for g in groups]
    flat = [i for g in groups for i in g]

    size = p.fov / 60 * 1.5  # cutout a bit larger than the grid so rotated footprints still fill it
    _stage(j, "download", f"Скачано 0 из {len(flat)} снимков", 0, len(flat))
    paths, failed = {}, 0
    with ThreadPoolExecutor(6) as ex:
        futs = {ex.submit(irsa.cutout, meta[i]["uri"], ra, dec, size): i for i in flat}
        for f in as_completed(futs):
            try:
                paths[futs[f]] = f.result()
            except Exception as e:
                failed += 1
                log.warning("cutout failed: %s", e)
            j.done += 1
            j.message = f"Скачано {j.done} из {j.total} снимков"

    if flat and failed == len(flat):
        raise ApiError(502, "upstream", f"IRSA недоступен: не удалось скачать ни одного снимка из {len(flat)}")

    _stage(j, "process", "Обработка", 0, len(groups))
    frames = []
    for g in groups:
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
        res["_partial"] = True

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
    d.mkdir(parents=True, exist_ok=True)
    vmin, vmax = prender.shared_limits([f.image for f in frames]) if frames else (0.0, 1.0)
    h = shape[0]
    for i, f in enumerate(frames):
        signed = p.mode == "diff"
        arr = diffs[i] if signed else f.image
        if signed:
            _atomic_png(prender.save_diff, d / f"f{i}.png", arr)
        else:
            _atomic_png(prender.save_gray, d / f"f{i}.png", arr, vmin, vmax, saturated=f.contam)
        _atomic_png(prender.save_mask, d / f"f{i}_flags.png", f.contam)  # bad and contaminated pixels; never SOURCE
        pred = None
        if r["motion"]:  # where the object should be on the frame date
            x, y = wcs.world_to_pixel(_sky(r, f.mjd))
            pred = prender.fits_to_png_xy(x, y, h)
        m = centroid.measure(arr, signed, saturated=f.contam)
        cands = [{"kind": c["kind"], "sign": c["sign"], "xy": prender.fits_to_png_xy(c["xy"][0], c["xy"][1], h),
                  "err_px": [round(e, 3) for e in c["err_px"]], "area_px": c["area"], "snr": round(c["snr"], 1),
                  "quality": c["quality"], "flags": c["flags"]}
                 for c in m["candidates"]]
        res["frames"].append({
            "index": i, "label": f"Эпоха {i + 1}", "mjd": round(f.mjd, 4), "date": f.date,
            "mjd_span": [round(f.mjd_min, 4), round(f.mjd_max, 4)], "det": f.detector,
            "wave_um": round(f.wave_um, 4) if np.isfinite(f.wave_um) else None, "n_exp": f.n_exp,
            "image_url": f"/cache/{j.key}/f{i}.png", "flags_url": f"/cache/{j.key}/f{i}_flags.png",
            "predicted_px": pred,
            # measured_px: the centre when the choice is unambiguous (status ok); null when nothing is found or
            # when the largest and the central object differ (both are listed, nothing is guessed). In diff mode
            # only the "+" residual (where the object is in this frame) gives the position.
            "measured_px": prender.fits_to_png_xy(m["xy"][0], m["xy"][1], h) if m["status"] == "ok" else None,
            "measured_candidates": cands,
            "measurement": {"status": m["status"], "reason": m["reason"], "quality": m["quality"],
                            "flags": m["flags"], "n_rejected_hot": m["n_rejected_hot"],
                            "n_low_snr": m["n_low_snr"],
                            "err_px": None if m["err_px"] is None else [round(e, 3) for e in m["err_px"]]},
            "span_days": round(f.span_days, 3), "bg_offsets": f.bg_offsets,
            "diff_noise_scale": scales[i] if scales else None,
            "obs_ids": f.obs_ids})
    if diffs is not None:
        _atomic_png(prender.save_rgb, d / "trail.png", trail(diffs), background=template, vmin=vmin, vmax=vmax)
        res["trail_url"] = f"/cache/{j.key}/trail.png"

    waves = [f.wave_um for f in frames if np.isfinite(f.wave_um)]
    if len(waves) > 1 and (max(waves) - min(waves)) / np.mean(waves) > WAVE_SPREAD_WARN:
        res["warnings"].append(f"λ в точке цели меняется между эпохами ({min(waves):.2f}-{max(waves):.2f} мкм): "
                               "у объектов с крутым спектром остатки в разностном кадре могут быть не движением")
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

    res = build_spectrum(center, lambda mjd: _sky(r, mjd), progress=progress)
    res["query"] = {"ra": ra, "dec": dec, "target": r["label"]}
    return res


def spectrum_view(p: Params):
    """Only target / ra+dec of p are used."""
    r = resolve.for_params(p)
    return _poll(cache.key_spectrum(*r["center"]), lambda j: _spectrum_job(j, r))


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
    return 200, res
