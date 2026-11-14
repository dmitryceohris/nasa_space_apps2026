"""Spectrum on the fly and light curve (spec sections 3.5, 3.6), built on src/Exposure + src/processing.

Exposures and cutouts come from this API's irsa module (same TAP lookup, same FITS cache as /api/view).
Photometry is done on the NATIVE pixels of every exposure (no reprojection, which would smear the flux).
Because of SPHEREx's linear variable filter each exposure sees the point at its own wavelength, read from the
WCS-WAVE extension, so the points of all 6 detectors together form the spectrum.
"""
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from astropy import units as u
from astropy.coordinates import SkyCoord

from Exposure import DETECTOR_BANDS, load_exposure
from processing.photometry import forced_photometry, light_curve

from . import irsa, resolve

log = logging.getLogger("spectrum")
CREDITS = "NASA SPHEREx via NASA/IPAC IRSA"
CUTOUT_DEG = 0.035  # ~20 px: aperture (2 px) + background annulus (4..7 px) with margin


def run(j, ra: float, dec: float, r: dict, per_epoch: int, stage) -> dict:
    """Background job body. r = resolve.for_params(...): if it carries proper motion, the aperture follows
    the object to its position at each exposure's date."""
    motion = r.get("motion")

    def coord_at(mjd: float) -> SkyCoord:
        pr, pd = resolve.propagate(motion, mjd) if motion else (ra, dec)
        return SkyCoord(pr * u.deg, pd * u.deg)

    stage(j, "query", "Поиск снимков всех 6 детекторов")
    chosen = []
    for det in DETECTOR_BANDS:
        for c in irsa.select(irsa.exposures(ra, dec, det), per_epoch):
            chosen += c

    stage(j, "download", f"Скачано 0 из {len(chosen)} снимков", 0, len(chosen))
    exps, failed = [], 0
    with ThreadPoolExecutor(6) as ex:
        futs = {ex.submit(irsa.cutout, e["uri"], ra, dec, CUTOUT_DEG): e for e in chosen}
        for f in as_completed(futs):
            try:
                exps.append(load_exposure(str(f.result())))
            except Exception as e:  # one lost exposure only removes one point of the spectrum
                failed += 1
                log.warning("spectrum exposure failed: %s", e)
            j.done += 1
            j.message = f"Скачано {j.done} из {j.total} снимков"

    stage(j, "process", "Фотометрия", 0, len(exps))
    pts = []
    for e in exps:
        p = forced_photometry(e, coord_at(e.mjd))
        if p is not None:
            pts.append(p)
        j.done += 1
    pts.sort(key=lambda p: p["wave_um"])

    warnings = []
    if failed:
        warnings.append(f"Не удалось получить {failed} из {len(chosen)} снимков")
    bad = sum(1 for p in pts if not p["ok"])
    if bad:
        warnings.append(f"{bad} точек с плохими пикселями в апертуре (часто насыщение яркой звезды): ok=false")
    return {"status": "ready", "query": {"ra": ra, "dec": dec, "target": r.get("label"), "per_epoch": per_epoch},
            "points": pts, "n_rejected": bad, "warnings": warnings, "credits": CREDITS,
            "bands_um": {str(k): [v[0], v[1]] for k, v in DETECTOR_BANDS.items()}}


def window(det: int | None, wave_min: float | None, wave_max: float | None) -> tuple[float, float]:
    """Narrow wavelength window for a light curve: explicit [wave_min, wave_max], else ±5% around the middle
    of detector `det`. A narrow window matters: with the LVF a wide one would mix variability with spectral shape."""
    if wave_min is not None and wave_max is not None:
        return float(wave_min), float(wave_max)
    lo, hi, _ = DETECTOR_BANDS[det or 2]
    mid = 0.5 * (lo + hi)
    return round(mid * 0.95, 4), round(mid * 1.05, 4)


def lightcurve(spec_result: dict, lo: float, hi: float) -> dict:
    pts = light_curve(spec_result["points"], lo, hi)
    return {"status": "ready", "query": spec_result.get("query"), "wave_window_um": [lo, hi],
            "points": [{k: p[k] for k in ("mjd", "flux_ujy", "err_ujy", "wave_um", "det", "ok", "flags")} for p in pts],
            "credits": CREDITS}
