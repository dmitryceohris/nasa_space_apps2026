"""Run: uvicorn api.main:app --port 8000   (or: bash start.sh)"""
import logging
import os
import threading
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response, RedirectResponse
from fastapi.staticfiles import StaticFiles

from . import cache, irsa, jobs, params, resolve
from .errors import ApiError, PUBLIC_INTERNAL_MESSAGE
from .limits import RequestLimiter
from runtime_storage import StorageFull, STORAGE_MESSAGE

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
VERSION = "0.1.0"
app = FastAPI(title="SPHEREx blink API", version=VERSION)


request_limiter = RequestLimiter(int(os.environ.get("API_RATE_LIMIT_PER_MINUTE", "120")))
_LIMITED_PATHS = {"/api/view", "/api/spectrum", "/api/lightcurve", "/api/resolve", "/api/export.gif"}


def _err(status, code, message, headers=None):
    return JSONResponse({"error": {"code": code, "message": message}}, status_code=status, headers=headers)


@app.middleware("http")
async def limit_api_requests(request: Request, call_next):
    if request.method != "OPTIONS" and request.url.path.rstrip("/") in _LIMITED_PATHS:
        peer = request.client.host if request.client else "unknown"
        retry = request_limiter.retry_after(peer)
        if retry is not None:
            return _err(429, "rate_limited", "Слишком много запросов. Повторите позже.",
                        {"Retry-After": str(retry)})
    return await call_next(request)


app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"],
                   allow_headers=["*"], expose_headers=["Retry-After"])


@app.exception_handler(ApiError)
async def _api_error(_: Request, e: ApiError):
    headers = {"Retry-After": "5"} if e.status == 503 else None
    return _err(e.status, e.code, PUBLIC_INTERNAL_MESSAGE if e.status == 500 else e.message, headers)


@app.exception_handler(StorageFull)
async def _storage_full(_: Request, e: StorageFull):
    logging.getLogger("api").warning("scratch capacity exhausted: %s", e)
    return _err(503, "storage_full", STORAGE_MESSAGE, {"Retry-After": "5"})


@app.exception_handler(RequestValidationError)
async def _validation_error(_: Request, e: RequestValidationError):
    bad = ", ".join(str(x["loc"][-1]) for x in e.errors())
    return _err(400, "bad_param", f"неверный или отсутствующий параметр: {bad}")


@app.exception_handler(StarletteHTTPException)
async def _http_error(_: Request, e: StarletteHTTPException):
    code = {404: "not_found", 405: "method_not_allowed"}.get(e.status_code, "http_error")
    return _err(e.status_code, code, PUBLIC_INTERNAL_MESSAGE if e.status_code >= 500 else str(e.detail), e.headers)


@app.exception_handler(Exception)
async def _unhandled(_: Request, e: Exception):
    logging.getLogger("api").exception("unhandled")
    return _err(500, "internal", PUBLIC_INTERNAL_MESSAGE, {"Access-Control-Allow-Origin": "*"})


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
         det: str | None = None, mode: str | None = None, max_per_epoch: str | None = None, snapshot: str | None = None):
    p = _p(target, ra, dec, fov, det, mode, max_per_epoch)
    if snapshot:
        return jobs.snapshot_result(p, snapshot)
    status, payload = jobs.view(p)
    return JSONResponse(payload, status_code=status)


@app.get("/api/export.gif")
def export_gif(target: str | None = None, ra: str | None = None, dec: str | None = None, fov: str | None = None,
               det: str | None = None, mode: str | None = None, max_per_epoch: str | None = None, snapshot: str | None = None):
    """GIF from already rendered frames. 409 if /api/view has not produced them yet."""
    p = _p(target, ra, dec, fov, det, mode, max_per_epoch)
    return Response(jobs.gif(p, snapshot) if snapshot else jobs.gif(p), media_type="image/gif")


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


_ROOT = Path(__file__).resolve().parent.parent


@app.get("/dev", include_in_schema=False)
def dev_redirect():
    return RedirectResponse("/dev/")


if os.environ.get("SPHEREX_DATABASE_URL"):
    @app.get("/cache/{key}/{filename}", include_in_schema=False)
    def stored_image(key: str, filename: str):
        from . import postgres
        data = postgres.read_image(key, filename)
        if data is None:
            raise ApiError(404, "not_found", "Кадр не найден")
        return Response(data, media_type="image/png")
else:
    app.mount("/cache", StaticFiles(directory=cache.CACHE), name="cache")
# "/" is the team's frontend (frontend/index.html, same origin as the API, so one link for the demo);
# the API developer's test page stays available at /dev/.
app.mount("/dev", StaticFiles(directory=_ROOT / "static", html=True), name="dev")
app.mount("/", StaticFiles(directory=_ROOT / "frontend", html=True), name="frontend")
