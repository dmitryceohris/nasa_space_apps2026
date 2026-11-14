"""SPHEREx Sky Viewer — бэкенд (FastAPI).

Запуск из корня репозитория:
    uv run uvicorn api.main:app --port 8000
Спецификация: docs/API.md
"""
from __future__ import annotations

import logging
import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import astropy.units as u  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.time import Time  # noqa: E402
from fastapi import FastAPI, Query, Request  # noqa: E402
from fastapi.middleware.cors import CORSMiddleware  # noqa: E402
from fastapi.responses import JSONResponse  # noqa: E402
from fastapi.staticfiles import StaticFiles  # noqa: E402

from Exposure import DETECTOR_BANDS  # noqa: E402
from processing import (Target, build_spectrum, build_view, find_preset,  # noqa: E402
                        lightcurve_from_spectrum, load_targets, resolve_simbad)

from .jobs import CACHE_DIR, JobManager, cache_key  # noqa: E402

VERSION = "0.2.0"
RESULT_VERSION = 2          # увеличить, если меняется формат результата — старый кэш станет невидим
CENTER_EPOCH = Time(2026.0, format="jyear")  # середина обзора SPHEREx: центр поля для движущихся целей
MAX_FOV = 20.0

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
log = logging.getLogger("spherex.api")

TARGETS = load_targets()
jobs = JobManager(workers=int(os.environ.get("JOB_WORKERS", "2")))


class ApiError(Exception):
    def __init__(self, status: int, code: str, message: str):
        self.status, self.code, self.message = status, code, message


# ------------------------------------------------------------------ параметры
def _parse_coord(ra: str, dec: str) -> SkyCoord:
    try:
        sexa = any(ch in ra for ch in "hms: ")
        if sexa:
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        else:
            c = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
    except Exception as e:  # noqa: BLE001
        raise ApiError(400, "bad_param", f"Не удалось разобрать координаты ra={ra!r} dec={dec!r}: {e}") from e
    return c


def _resolve_target(name: str) -> Target:
    t = find_preset(name, TARGETS)
    if t is not None:
        return t
    try:
        t = resolve_simbad(name)
    except Exception as e:  # noqa: BLE001
        raise ApiError(502, "upstream", f"SIMBAD недоступен: {e}") from e
    if t is None:
        raise ApiError(404, "not_found", f"Объект {name!r} не найден в SIMBAD")
    return t


def resolve_params(target: str | None, ra: str | None, dec: str | None) -> tuple[SkyCoord, Target | None]:
    """Центр поля и (если известна) цель с собственным движением."""
    tgt = _resolve_target(target) if target else None
    if ra is not None and dec is not None:
        center = _parse_coord(ra, dec)
        if not (0 <= center.ra.deg <= 360 and -90 <= center.dec.deg <= 90):
            raise ApiError(400, "bad_param", "ra должен быть в 0..360, dec в -90..90")
        return center, tgt
    if ra is not None or dec is not None:
        raise ApiError(400, "bad_param", "Нужны оба параметра: ra и dec")
    if tgt is None:
        raise ApiError(400, "bad_param", "Укажите target или пару ra+dec")
    return tgt.coord_at(CENTER_EPOCH), tgt


def _check_range(name, value, lo, hi):
    if not (lo <= value <= hi):
        raise ApiError(400, "bad_param", f"{name} должен быть в диапазоне {lo}..{hi}")


def _reply(status: int, body: dict) -> JSONResponse:
    return JSONResponse(body, status_code=status)


# ------------------------------------------------------------------ приложение
@asynccontextmanager
async def lifespan(_app: FastAPI):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if os.environ.get("PREWARM", "1") != "0":
        # прогрев: на демо ничего не должно качаться вживую
        for t in TARGETS.values():
            try:
                _start_view(t.coord_at(CENTER_EPOCH), t, t.fov, t.det, 4)
            except Exception as e:  # noqa: BLE001
                log.warning("prewarm %s: %s", t.id, e)
    yield


app = FastAPI(title="SPHEREx Sky Viewer", version=VERSION, lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.exception_handler(ApiError)
async def _api_error(_req: Request, e: ApiError):
    return JSONResponse({"error": {"code": e.code, "message": e.message}}, status_code=e.status)


@app.exception_handler(Exception)
async def _internal(_req: Request, e: Exception):
    log.exception("internal error")
    return JSONResponse({"error": {"code": "internal", "message": f"{type(e).__name__}: {e}"}}, status_code=500)


@app.get("/api/health")
def health():
    return {"status": "ok", "version": VERSION}


@app.get("/api/targets")
def targets():
    return [t.to_json() for t in TARGETS.values()]


@app.get("/api/resolve")
def resolve(target: str):
    t = _resolve_target(target)
    now = t.coord_at(Time.now())
    return {**t.to_json(), "ra_now": now.ra.deg, "dec_now": now.dec.deg}


def _start_view(center: SkyCoord, tgt: Target | None, fov: float, det: int, max_per_epoch: int):
    key = cache_key(v=RESULT_VERSION, ra=round(center.ra.deg, 5), dec=round(center.dec.deg, 5),
                    fov=fov, det=det, n=max_per_epoch, t=(tgt.name if tgt else None))
    prefix = f"/cache/view/{key}"

    def run(out_dir: Path, progress):
        res = build_view(center, fov, det, max_per_epoch, out_dir, prefix, target=tgt, progress=progress)
        res["query"] = {"ra": round(center.ra.deg, 6), "dec": round(center.dec.deg, 6), "fov": fov,
                        "det": det, "max_per_epoch": max_per_epoch,
                        "target": tgt.to_json() if tgt else None,
                        "band_um": list(DETECTOR_BANDS[det][:2])}
        return res

    return jobs.get_or_start("view", key, run)


@app.get("/api/view")
def view(target: str | None = None, ra: str | None = None, dec: str | None = None,
         fov: float = 3.0, det: int = 3, max_per_epoch: int = Query(4, ge=1, le=20)):
    _check_range("fov", fov, 1, MAX_FOV)
    _check_range("det", det, 1, 6)
    center, tgt = resolve_params(target, ra, dec)
    return _reply(*_start_view(center, tgt, fov, det, max_per_epoch))


def _start_spectrum(center: SkyCoord, tgt: Target | None, per_epoch: int):
    key = cache_key(v=RESULT_VERSION, ra=round(center.ra.deg, 4), dec=round(center.dec.deg, 4),
                    n=per_epoch, t=(tgt.name if tgt else None))
    coord_at = tgt.coord_at if tgt else (lambda _mjd: center)

    def run(_out_dir: Path, progress):
        res = build_spectrum(center, coord_at, per_epoch=per_epoch, progress=progress)
        res["query"] = {"ra": round(center.ra.deg, 6), "dec": round(center.dec.deg, 6),
                        "target": tgt.name if tgt else None}
        return res

    return jobs.get_or_start("spectrum", key, run)


@app.get("/api/spectrum")
def spectrum_ep(target: str | None = None, ra: str | None = None, dec: str | None = None,
                per_epoch: int = Query(2, ge=1, le=10)):
    center, tgt = resolve_params(target, ra, dec)
    # клик по картинке передаёт ra/dec точки; цель с собственным движением учитываем,
    # только если точка совпадает с ней (иначе апертура «поедет» за чужим объектом)
    if tgt is not None and ra is not None and center.separation(tgt.coord_at(CENTER_EPOCH)).arcsec > 12:
        tgt = None
    return _reply(*_start_spectrum(center, tgt, per_epoch))


@app.get("/api/lightcurve")
def lightcurve(target: str | None = None, ra: str | None = None, dec: str | None = None,
               det: int = 2, wave_min: float | None = None, wave_max: float | None = None,
               per_epoch: int = Query(2, ge=1, le=10)):
    _check_range("det", det, 1, 6)
    lo, hi, _ = DETECTOR_BANDS[det]
    mid = 0.5 * (lo + hi)
    wmin = wave_min if wave_min is not None else mid * 0.95
    wmax = wave_max if wave_max is not None else mid * 1.05
    center, tgt = resolve_params(target, ra, dec)
    if tgt is not None and ra is not None and center.separation(tgt.coord_at(CENTER_EPOCH)).arcsec > 12:
        tgt = None
    status, body = _start_spectrum(center, tgt, per_epoch)
    if status != 200 or body.get("status") != "ready":
        return _reply(status, body)
    return lightcurve_from_spectrum(body, wmin, wmax)


# ------------------------------------------------------------------ статика
CACHE_DIR.mkdir(parents=True, exist_ok=True)
app.mount("/cache", StaticFiles(directory=CACHE_DIR), name="cache")
app.mount("/", StaticFiles(directory=ROOT / "frontend", html=True), name="frontend")
