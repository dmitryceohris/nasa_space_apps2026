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
    name = next((x["name"] for x in resolve.PRESETS if x["id"] == target.lower()), target)
    return resolve.resolve(name)


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


# Not implemented: left empty on purpose (nothing is better than a wrong number). Contract: spec sections 3.5, 3.6.
@app.get("/api/spectrum")
def spectrum(ra: str | None = None, dec: str | None = None, target: str | None = None):
    # TODO: photometry at the point over all exposures of all 6 detectors -> {"points":[{wave_um,bw_um,flux_ujy,
    # err_ujy,mjd,det,flags}], "n_rejected"}. Same 202/200 job mechanism as /api/view.
    raise ApiError(501, "not_implemented", "спектр ещё не реализован")


@app.get("/api/lightcurve")
def lightcurve(ra: str | None = None, dec: str | None = None, wave_min: float | None = None,
               wave_max: float | None = None):
    # TODO: flux vs time in a narrow wavelength window -> {"wave_window_um":[lo,hi], "points":[{mjd,flux_ujy,err_ujy}]}
    raise ApiError(501, "not_implemented", "кривая блеска ещё не реализована")


@app.on_event("startup")
def _startup():
    if os.environ.get("PREWARM", "1") == "1":
        threading.Thread(target=jobs.prewarm, daemon=True).start()


app.mount("/cache", StaticFiles(directory=cache.CACHE), name="cache")
app.mount("/", StaticFiles(directory=Path(__file__).resolve().parent.parent / "static", html=True), name="static")
