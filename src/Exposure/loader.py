"""Загрузка FITS-файла SPHEREx (полного или вырезки) в Exposure."""
from __future__ import annotations

import os
import threading
import warnings

import numpy as np
from astropy import log
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

from .exposure import Exposure, parse_flag_bits
from runtime_storage import current_workspace

# Ключи времени. В реальных файлах SPHEREx (проверено на DR1) есть MJD-BEG/MJD-AVG/MJD-END —
# берём MJD-AVG (середина экспозиции). Остальные — запасные варианты.
_MJD_KEYS = ("MJD-AVG", "MJD-MID", "MJDMID", "MJD-OBS", "MJD-BEG", "MJD_OBS", "MJD")
_DATE_KEYS = ("DATE-AVG", "DATE-MID", "DATE-OBS", "DATE-BEG")
# ВНИМАНИЕ: "DATE" — это дата записи файла, а не наблюдения. Её не используем.


def _find_hdu(hdul: fits.HDUList, *names: str):
    for n in names:
        if n in hdul:
            return hdul[n]
    return None


def _read_time_info(headers, mjd_hint: float | None):
    # Midpoint timestamps take priority even if they occur in the primary header.
    for keys, fmt in ((_MJD_KEYS[:1], "mjd"), (_DATE_KEYS[:1], None),
                      (_MJD_KEYS[1:], "mjd"), (_DATE_KEYS[1:], None)):
        for key in keys:
            for h in headers:
                if key in h:
                    declared_scale = next((str(header["TIMESYS"]).strip().lower()
                                           for header in [h] + [other for other in headers if other is not h]
                                           if str(header.get("TIMESYS", "")).strip()), None)
                    assumptions = [] if declared_scale else ["timesys_assumed_utc"]
                    if "BEG" in key or "OBS" in key or key == "MJD":
                        assumptions.append("timestamp_not_verified_midpoint")
                    value = float(h[key]) if fmt else h[key]
                    return float(Time(value, format=fmt, scale=declared_scale or "utc").utc.mjd), key, assumptions
    if mjd_hint is not None:  # время из ответа архива (t_min/t_max)
        return float(mjd_hint), "archive_hint", ["archive_time_may_be_lower_bound"]
    raise ValueError("Не найдено время наблюдения ни в заголовке, ни в mjd_hint")


def _read_mjd(headers, mjd_hint: float | None) -> float:
    return _read_time_info(headers, mjd_hint)[0]


def _read_detector(headers, source: str) -> int:
    for h in headers:
        if "DETECTOR" in h:
            return int(h["DETECTOR"])
    # запасной вариант: имя файла ..._2D3_spx_l2b... -> 3
    import re
    m = re.search(r"D([1-6])_spx", source)
    if m:
        return int(m.group(1))
    raise ValueError("Не удалось определить номер детектора")


def exposure_from_hdul(hdul: fits.HDUList, source: str = "", *,
                       mjd_hint: float | None = None, obs_id_hint: str = "") -> Exposure:
    img_hdu = _find_hdu(hdul, "IMAGE")
    if img_hdu is None:
        raise ValueError(f"{source}: нет расширения IMAGE — это не спектральный снимок SPHEREx")
    flags_hdu = _find_hdu(hdul, "FLAGS")
    var_hdu = _find_hdu(hdul, "VARIANCE")
    zodi_hdu = _find_hdu(hdul, "ZODI")

    image = np.asarray(img_hdu.data, dtype=np.float32)
    if image.ndim != 2:
        raise ValueError(f"{source}: IMAGE must be two-dimensional")
    if var_hdu is None:  # without variance every pixel would be marked bad and the frame silently empty
        raise ValueError(f"{source}: нет расширения VARIANCE — без дисперсии нельзя маскировать и взвешивать пиксели")
    variance = np.asarray(var_hdu.data, dtype=np.float32)
    flags = (np.asarray(flags_hdu.data, dtype=np.int32) if flags_hdu is not None
             else np.zeros(image.shape, dtype=np.int32))
    zodi = np.asarray(zodi_hdu.data, dtype=np.float32) if zodi_hdu is not None else None
    if any(a is not None and a.shape != image.shape for a in (variance, flags, zodi)):
        raise ValueError(f"{source}: IMAGE/FLAGS/VARIANCE/ZODI shape mismatch")

    hdr = img_hdu.header
    headers = [hdr, hdul[0].header]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(hdr)  # пространственный WCS (TAN-SIP)
        old_level = log.level
        log.setLevel("WARNING")  # прячем INFO про «Inconsistent SIP» — известная особенность key='W'
        try:
            # Спектральный WCS: альтернативный ключ 'W' + таблица WCS-WAVE (как в туториалах IRSA)
            wave_wcs = WCS(header=hdr, fobj=hdul, key="W")
            wave_wcs.sip = None
        except Exception as e:  # noqa: BLE001 — без λ всё равно можно показать картинку
            warnings.warn(f"{source}: не удалось прочитать WCS-WAVE ({e}); λ будет приблизительной")
            wave_wcs = None
        finally:
            log.setLevel(old_level)

    mjd, time_source, time_assumptions = _read_time_info(headers, mjd_hint)
    epsf_hdu = _find_hdu(hdul, "EPSF")
    exp = Exposure(
        image=image, variance=variance, flags=flags, zodi=zodi,
        wcs=wcs, wave_wcs=wave_wcs,
        detector=_read_detector(headers, source),
        mjd=mjd, time_source=time_source, time_assumptions=time_assumptions,
        obs_id=str(hdr.get("OBSID", hdul[0].header.get("OBSID", obs_id_hint))),
        flag_bits=parse_flag_bits(flags_hdu.header) if flags_hdu is not None else {},
        source=source, header=hdr,
        flags_available=flags_hdu is not None,
        epsf_data=epsf_hdu.data.copy() if epsf_hdu is not None else None,
        epsf_header=epsf_hdu.header.copy() if epsf_hdu is not None else None,
    )
    from processing.observations import product_record
    from .calibration import attach_registered_calibration
    exp.provenance = product_record(source, obs_id=exp.obs_id, detector=exp.detector)
    exp.provenance["pipeline_version"] = str(hdul[0].header.get("VERSION", "unknown"))
    return attach_registered_calibration(exp, exp.provenance["pipeline_version"])


# Общий на весь процесс лимит одновременных загрузок с IRSA. При ~12 параллельных запросах
# сервис начинает массово отвечать 503 (проверено), поэтому все потоки и задачи делят этот семафор.
_IRSA_SLOTS = threading.BoundedSemaphore(int(os.environ.get("IRSA_MAX_CONCURRENT", "4")))


def _download(url: str, cache: bool, timeout: float, retries: int = 5) -> str:
    """Скачивание с повторами и экспоненциальной паузой: IRSA периодически отвечает 503."""
    import random
    import time
    from urllib.error import HTTPError, URLError

    from astropy.utils.data import download_file, is_url_in_cache
    if cache and is_url_in_cache(url):
        return download_file(url, cache=True)  # из кэша — без очереди к IRSA
    for attempt in range(retries):
        try:
            with _IRSA_SLOTS:
                return download_file(url, cache=cache, timeout=timeout)
        except (HTTPError, URLError, TimeoutError, ConnectionError) as e:
            code = getattr(e, "code", None)
            if attempt == retries - 1 or (code is not None and code not in (429, 500, 502, 503, 504)):
                raise
            time.sleep(2 ** (attempt + 1) + random.uniform(0, 1.5))
    raise RuntimeError("unreachable")


def load_exposure(path_or_url: str, *, mjd_hint: float | None = None,
                  obs_id_hint: str = "", cache: bool = True, timeout: float = 60,
                  product_url: str | None = None, workspace=None) -> Exposure:
    """Открывает локальный файл или URL (URL кэшируется на диске средствами astropy)."""
    if workspace is not None:
        # ContextVars do not cross the download worker's thread boundary automatically.
        with workspace.activate():
            return load_exposure(path_or_url, mjd_hint=mjd_hint, obs_id_hint=obs_id_hint,
                                 cache=cache, timeout=timeout, product_url=product_url)
    workspace = current_workspace()
    src = path_or_url
    owned = False
    if path_or_url.startswith(("http://", "https://")):
        if workspace is not None:
            with _IRSA_SLOTS:
                src = workspace.download(path_or_url, timeout=timeout)
            owned = True
        else:
            src = _download(path_or_url, cache, timeout)
    try:
        with fits.open(src, memmap=False) as hdul:
            exp = exposure_from_hdul(hdul, source=product_url or path_or_url,
                                     mjd_hint=mjd_hint, obs_id_hint=obs_id_hint)
        import hashlib
        from pathlib import Path
        exp.provenance["file_sha256"] = hashlib.sha256(Path(src).read_bytes()).hexdigest()
        return exp
    finally:
        if owned:
            src.unlink(missing_ok=True)
