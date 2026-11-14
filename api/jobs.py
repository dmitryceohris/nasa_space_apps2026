"""Background jobs, one per cache key. /api/view starts a job and the front polls the same URL until it is ready."""
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
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from reproject import reproject_interp

from . import cache, centroid, irsa, render, resolve
from .errors import ApiError
from .params import Params

log = logging.getLogger("jobs")
CREDITS = "NASA SPHEREx via NASA/IPAC IRSA"
POOL = ThreadPoolExecutor(max_workers=4)
JOBS = {}
LOCK = threading.Lock()
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


def view(p: Params):
    """-> (http_status, payload): 200 cached/ready/error, 202 running."""
    if p.mode == "trail":
        raise ApiError(501, "not_implemented", "режим trail ещё не реализован")
    r = resolve.for_params(p)
    ra, dec = r["center"]
    k = cache.key(ra, dec, p.fov, p.det, p.mode, p.max_per_epoch)
    res = cache.read(k)
    if res:
        return 200, res
    with LOCK:
        j = JOBS.get(k)
        if j is None:
            j = JOBS[k] = Job(k, uuid.uuid4().hex[:8])
            POOL.submit(_run, j, p, r)
        if j.status == "error":  # drop it so the next request retries
            JOBS.pop(k)
            return 200, {"status": "error", "message": j.message}
        if j.status == "ready":
            JOBS.pop(k)
            return 200, j.result
        return 202, {"status": "running", "job_id": j.id, "stage": j.stage, "done": j.done, "total": j.total,
                     "message": j.message}


def cached_result(p: Params):
    r = resolve.for_params(p)
    ra, dec = r["center"]
    return cache.read(cache.key(ra, dec, p.fov, p.det, p.mode, p.max_per_epoch))


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


def _run(j, p, r):
    try:
        res = _pipeline(j, p, r)
        if res["frames"]:
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


def _pipeline(j, p, r):
    ra, dec = r["center"]
    w, shape = render.grid(ra, dec, p.fov)
    res = {"status": "ready",
           "query": {"ra": ra, "dec": dec, "fov": p.fov, "det": p.det, "mode": p.mode, "target": r["label"]},
           "grid": {"width_px": shape[1], "height_px": shape[0], "pixel_scale_arcsec": render.PIX_ARCSEC,
                    "wcs": render.wcs_dict(w)},
           "frames": [], "trail_url": None, "motion": None, "warnings": [], "credits": CREDITS}

    _stage(j, "query", "Поиск снимков")
    exps = irsa.exposures(ra, dec, p.det)
    if not exps:
        res["warnings"].append("SPHEREx ещё не наблюдал эту область")
        return res
    clusters = irsa.select(exps, p.max_per_epoch)
    if p.mode == "single":  # one exposure per epoch: download only those
        clusters = [c[:1] for c in clusters]
    flat = [e for c in clusters for e in c]

    size = p.fov / 60 * 1.5  # cutout a bit larger than the grid so rotated footprints still fill it
    _stage(j, "download", f"Скачано 0 из {len(flat)} снимков", 0, len(flat))
    paths, failed = {}, 0
    with ThreadPoolExecutor(6) as ex:
        futs = {ex.submit(irsa.cutout, e["uri"], ra, dec, size): e["uri"] for e in flat}
        for f in as_completed(futs):
            try:
                paths[futs[f]] = f.result()
            except Exception as e:
                failed += 1
                log.warning("cutout failed: %s", e)
            j.done += 1
            j.message = f"Скачано {j.done} из {j.total} снимков"

    _stage(j, "process", "Обработка", 0, len(clusters))
    epochs = []
    for c in clusters:
        imgs, flags, used = [], np.zeros(shape, bool), []
        for e in c:
            if e["uri"] not in paths:
                continue
            try:
                with fits.open(paths[e["uri"]]) as h:
                    im, _ = reproject_interp((h["IMAGE"].data, WCS(h["IMAGE"].header, fobj=h)), w, shape)
                    fl, _ = reproject_interp((h["FLAGS"].data, WCS(h["FLAGS"].header, fobj=h)), w, shape,
                                             order="nearest-neighbor")
            except Exception as ex_:
                failed += 1
                log.warning("reproject failed: %s", ex_)
                continue
            imgs.append(im)
            flags |= np.nan_to_num(fl) != 0
            used.append(e)
        if used:
            epochs.append({"img": np.nanmean(imgs, axis=0), "flags": flags,
                           "mjd": float(np.mean([e["mjd"] for e in used])),
                           "wave": float(np.mean([sum(e["lam"]) / 2 for e in used])),
                           "n": len(used), "obs": [o for o in (irsa.obs_id(e["uri"]) for e in used) if o]})
        j.done += 1
    if failed:
        res["warnings"].append(f"Не удалось получить или обработать {failed} из {len(flat)} снимков")
    if p.mode == "diff":
        if len(epochs) < 2:
            res["warnings"].append("diff требует минимум 2 эпохи")
            epochs = []
        else:  # this epoch minus the median over all epochs
            tmpl = np.nanmedian([e["img"] for e in epochs], axis=0)
            for e in epochs:
                e["img"] = e["img"] - tmpl

    _stage(j, "render", "Рисую кадры", 0, len(epochs))
    signed = p.mode == "diff"
    lo, hi = render.limits([e["img"] for e in epochs], signed)
    d = cache.CACHE / j.key
    d.mkdir(exist_ok=True)
    for i, e in enumerate(epochs):
        render.save_png(e["img"], d / f"f{i}.png", lo, hi, signed)
        render.save_flags(e["flags"], d / f"f{i}_flags.png")
        pred = None
        if r["motion"] is not None:  # where the object should be on this date (proper motion + parallax)
            pr, pd = resolve.propagate(r["motion"], e["mjd"])
            x, y = w.world_to_pixel(SkyCoord(pr * u.deg, pd * u.deg))
            pred = render.png_xy(x, y, shape[0])
        cands = [{"kind": c["kind"], "xy": render.png_xy(c["xy"][0], c["xy"][1], shape[0]), "area_px": c["area"]}
                 for c in centroid.find(e["img"], signed)]
        res["frames"].append({
            "index": i, "label": f"Эпоха {i + 1}", "mjd": e["mjd"], "date": Time(e["mjd"], format="mjd").iso[:10],
            "det": p.det, "wave_um": e["wave"], "n_exp": e["n"],
            "image_url": f"/cache/{j.key}/f{i}.png", "flags_url": f"/cache/{j.key}/f{i}_flags.png",
            # measured_px: the centre when there is exactly one candidate; null when none or when the largest and
            # the central object differ (both are then listed in measured_candidates instead of guessing one).
            "predicted_px": pred, "measured_px": cands[0]["xy"] if len(cands) == 1 else None,
            "measured_candidates": cands, "obs_ids": e["obs"]})
    return res
