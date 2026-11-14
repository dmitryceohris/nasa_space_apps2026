"""Run: uvicorn api.main:app --port 8000   (or: bash start.sh)"""
import logging
import os
import threading
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from . import cache, irsa, jobs, params, resolve
from .errors import ApiError

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
VERSION = "0.1.0"
app = FastAPI(title="SPHEREx blink API", version=VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


def _err(status, code, message):
    return JSONResponse({"error": {"code": code, "message": message}}, status_code=status)


@app.exception_handler(ApiError)
async def _api_error(_: Request, e: ApiError):
    return _err(e.status, e.code, e.message)


@app.exception_handler(Exception)
async def _unhandled(_: Request, e: Exception):
    logging.getLogger("api").exception("unhandled")
    return _err(500, "internal", str(e))


def _p(target, ra, dec, fov, det, mode, max_per_epoch):
    return params.parse(target, ra, dec, fov, det, mode, max_per_epoch)


@app.get("/api/health")
def health():
    return {"status": "ok", "version": VERSION}


@app.get("/api/targets")
def targets():
    return resolve.PRESETS


@app.get("/api/sources")
def sources():
    """Local exposure tables currently loaded (folder tables/, or SPHEREX_TABLES)."""
    return irsa.sources()


@app.get("/api/resolve")
def resolve_ep(target: str):
    return resolve.resolve(target)


@app.get("/api/view")
def view(target: str | None = None, ra: str | None = None, dec: str | None = None, fov: str | None = None,
         det: str | None = None, mode: str | None = None, max_per_epoch: str | None = None):
    status, payload = jobs.view(_p(target, ra, dec, fov, det, mode, max_per_epoch))
    return JSONResponse(payload, status_code=status)


@app.get("/api/export.gif")
def export_gif(target: str | None = None, ra: str | None = None, dec: str | None = None, fov: str | None = None,
               det: str | None = None, mode: str | None = None, max_per_epoch: str | None = None):
    """GIF from already rendered frames. 409 if /api/view has not produced them yet."""
    return Response(jobs.gif(_p(target, ra, dec, fov, det, mode, max_per_epoch)), media_type="image/gif")


@app.get("/api/spectrum")
def spectrum(target: str | None = None, ra: str | None = None, dec: str | None = None):
    """Photometry at the target over all exposures of all 6 detectors (src/processing.build_spectrum).
    Same 202/200 job mechanism as /api/view."""
    status, payload = jobs.spectrum_view(params.parse(target, ra, dec))
    return JSONResponse(payload, status_code=status)


@app.get("/api/lightcurve")
def lightcurve(target: str | None = None, ra: str | None = None, dec: str | None = None,
               wave_min: str | None = None, wave_max: str | None = None, det: str | None = None):
    """Flux over time in a narrow wavelength window of the (cached) spectrum.
    Window defaults to +-5% around the middle of detector `det` (default 3)."""
    p = params.parse(target, ra, dec, det=det)
    lo, hi = jobs.default_window(p.det)
    try:
        lo = float(wave_min) if wave_min not in (None, "") else lo
        hi = float(wave_max) if wave_max not in (None, "") else hi
    except ValueError:
        raise ApiError(400, "bad_param", "wave_min и wave_max должны быть числами (мкм)")
    if not 0 < lo < hi:
        raise ApiError(400, "bad_param", "нужно 0 < wave_min < wave_max")
    status, payload = jobs.lightcurve_view(p, lo, hi)
    return JSONResponse(payload, status_code=status)


@app.on_event("startup")
def _startup():
    if os.environ.get("PREWARM", "1") == "1":
        threading.Thread(target=jobs.prewarm, daemon=True).start()


app.mount("/cache", StaticFiles(directory=cache.CACHE), name="cache")
app.mount("/", StaticFiles(directory=Path(__file__).resolve().parent.parent / "static", html=True), name="static")
