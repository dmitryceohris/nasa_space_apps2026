"""Background jobs, one per cache key. The front polls the same URL until the job is ready (202 -> 200).

Image processing is NOT done here: the frames come from src/processing (coadd with zodiacal-light subtraction,
bad-pixel masks, wavelength at the target, difference, trail) and the spectrum from build_spectrum. This module
only finds the exposures (local table or IRSA TAP), downloads the cutouts, calls those functions, writes the
PNG files and adds the object centres.
"""
import io
import logging
import threading
import uuid
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from Exposure import DETECTOR_BANDS, load_exposure
from processing import render as prender
from processing.epochs import select_per_epoch
from processing.frames import coadd, difference, grid_header_json, make_grid, trail
from processing.pipeline import CREDITS, build_spectrum, lightcurve_from_spectrum

from . import cache, centroid, irsa, render, resolve
from .errors import ApiError
from .params import Params

log = logging.getLogger("jobs")
POOL = ThreadPoolExecutor(max_workers=4)
JOBS = {}
LOCK = threading.Lock()
WAVE_SPREAD_WARN = 0.05  # same threshold as src/processing/pipeline.py: lambda change between epochs
warnings.filterwarnings("ignore")


@dataclass
class Job:
    key: str
    id: str
    status: str = "running"  # running | ready (empty result, not cached) | error
    stage: str = "query"  # query -> download -> process -> render
    done: int = 0
    total: int = 0
    message: str = ""
    result: dict | None = None


def _poll(k: str, runner):
    """-> (http_status, payload): 200 cached/ready/error, 202 running. Starts runner(job) if nothing runs for k."""
    res = cache.read(k)
    if res:
        return 200, res
    with LOCK:
        j = JOBS.get(k)
        if j is None:
            j = JOBS[k] = Job(k, uuid.uuid4().hex[:8])
            POOL.submit(_run, j, runner)
        if j.status == "error":  # drop it so the next request retries
            JOBS.pop(k)
            return 200, {"status": "error", "message": j.message}
        if j.status == "ready":
            JOBS.pop(k)
            return 200, j.result
        return 202, {"status": "running", "job_id": j.id, "stage": j.stage, "done": j.done, "total": j.total,
                     "message": j.message}


def _run(j, runner):
    try:
        res = runner(j)
        empty = (res.get("frames") is not None and not res["frames"]) or (
            res.get("points") is not None and not res["points"])
        if not empty:
            cache.write(j.key, res)
            with LOCK:
                JOBS.pop(j.key, None)
        else:  # nothing to show: do not cache, new SPHEREx releases may fill it later
            j.result, j.status = res, "ready"
    except ApiError as e:
        j.status, j.message = "error", e.message
    except Exception as e:
        log.exception("job %s failed", j.key)
        j.status, j.message = "error", f"внутренняя ошибка: {e}"


def _stage(j, stage, msg="", done=0, total=0):
    j.stage, j.message, j.done, j.total = stage, msg, done, total


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
           "warnings": [], "credits": CREDITS}

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

    diffs = template = None
    if p.mode in ("diff", "trail"):
        if len(frames) < 2:
            res["warnings"].append(f"{p.mode} требует минимум 2 эпохи")
            if p.mode == "diff":  # nothing meaningful to show
                frames = []
        else:
            diffs, template = difference(frames)

    _stage(j, "render", "Рисую кадры", 0, len(frames))
    d = cache.CACHE / j.key
    d.mkdir(exist_ok=True)
    vmin, vmax = prender.shared_limits([f.image for f in frames]) if frames else (0.0, 1.0)
    h = shape[0]
    for i, f in enumerate(frames):
        signed = p.mode == "diff"
        arr = diffs[i] if signed else f.image
        if signed:
            prender.save_diff(arr, d / f"f{i}.png")
        else:
            prender.save_gray(arr, d / f"f{i}.png", vmin, vmax, saturated=f.contam)
        prender.save_mask(f.contam, d / f"f{i}_flags.png")  # bad and contaminated pixels; never SOURCE
        pred = None
        if r["motion"]:  # where the object should be on the frame date
            x, y = wcs.world_to_pixel(_sky(r, f.mjd))
            pred = prender.fits_to_png_xy(x, y, h)
        cands = [{"kind": c["kind"], "xy": prender.fits_to_png_xy(c["xy"][0], c["xy"][1], h), "area_px": c["area"]}
                 for c in centroid.find(arr, signed)]
        res["frames"].append({
            "index": i, "label": f"Эпоха {i + 1}", "mjd": round(f.mjd, 4), "date": f.date,
            "mjd_span": [round(f.mjd_min, 4), round(f.mjd_max, 4)], "det": f.detector,
            "wave_um": round(f.wave_um, 4) if np.isfinite(f.wave_um) else None, "n_exp": f.n_exp,
            "image_url": f"/cache/{j.key}/f{i}.png", "flags_url": f"/cache/{j.key}/f{i}_flags.png",
            "predicted_px": pred,
            # measured_px: the centre when there is exactly one candidate; null when none or when the largest
            # and the central object differ (both are then listed in measured_candidates, nothing is guessed).
            "measured_px": cands[0]["xy"] if len(cands) == 1 else None,
            "measured_candidates": cands, "obs_ids": f.obs_ids})
    if diffs is not None:
        prender.save_rgb(trail(diffs), d / "trail.png", background=template, vmin=vmin, vmax=vmax)
        res["trail_url"] = f"/cache/{j.key}/trail.png"

    waves = [f.wave_um for f in frames if np.isfinite(f.wave_um)]
    if len(waves) > 1 and (max(waves) - min(waves)) / np.mean(waves) > WAVE_SPREAD_WARN:
        res["warnings"].append(f"λ в точке цели меняется между эпохами ({min(waves):.2f}-{max(waves):.2f} мкм): "
                               "у объектов с крутым спектром остатки в разностном кадре могут быть не движением")
    c = (h // 2, shape[1] // 2)
    if any(f.contam[c[0] - 1:c[0] + 2, c[1] - 1:c[1] + 2].any() for f in frames):
        res["warnings"].append("В центре поля есть насыщенные или загрязнённые пиксели (слой флагов)")
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
