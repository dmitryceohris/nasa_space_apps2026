"""Общая сетка неба, сложение экспозиций по эпохам, разностные кадры и «след» движения."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.wcs import WCS

from Exposure import PIXEL_SCALE_ARCSEC, Exposure


# ---------------------------------------------------------------- сетка
def make_grid(center: SkyCoord, fov_arcmin: float,
              pix_arcsec: float = PIXEL_SCALE_ARCSEC) -> tuple[WCS, tuple[int, int]]:
    """Квадратная TAN-сетка: север вверх, восток влево, пиксель = родной пиксель SPHEREx.

    Мельче пиксель не делаем: это создало бы фальшивое разрешение.
    """
    n = max(5, int(round(fov_arcmin * 60 / pix_arcsec)))
    if n % 2 == 0:
        n += 1  # нечётный размер — центр ровно в пикселе
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [center.ra.deg, center.dec.deg]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]  # FITS 1-based
    w.wcs.cdelt = [-pix_arcsec / 3600, pix_arcsec / 3600]
    return w, (n, n)


def grid_header_json(wcs: WCS, shape) -> dict:
    return {"width_px": shape[1], "height_px": shape[0], "pixel_scale_arcsec": abs(wcs.wcs.cdelt[1]) * 3600,
            "wcs": {"CTYPE1": wcs.wcs.ctype[0], "CTYPE2": wcs.wcs.ctype[1],
                    "CRVAL1": float(wcs.wcs.crval[0]), "CRVAL2": float(wcs.wcs.crval[1]),
                    "CRPIX1": float(wcs.wcs.crpix[0]), "CRPIX2": float(wcs.wcs.crpix[1]),
                    "CDELT1": float(wcs.wcs.cdelt[0]), "CDELT2": float(wcs.wcs.cdelt[1])}}


SPIKE_SNR = 20.0  # a native pixel this many sigma above the sky ...
SPIKE_NEIGHBOUR_RATIO = 0.05  # ... whose brightest neighbour is below this fraction of it is a lone spike
# Real SPHEREx stars (1028 maxima with S/N > 30 in 160 cutouts, detectors 1-6) never had the brightest
# neighbour below 0.084 of the peak; a pixel-sized cosmic ray or hot pixel has neighbours at the noise level.


def lone_spike_mask(exp: Exposure, sub: np.ndarray | None = None, bad: np.ndarray | None = None) -> np.ndarray:
    """Unflagged single-pixel outliers (cosmic rays, hot pixels missed by the pipeline) in NATIVE pixels.

    Must be applied before interpolation: bilinear reprojection spreads one pixel over 2x2 grid pixels, after
    which it looks like a compact star (one such pixel became a S/N 279 "object" in an empty field).
    Needs at least 6 valid neighbours, so pixels next to masked (e.g. saturated) regions are never removed.
    sub / bad: exp.background_subtracted() / exp.bad_mask() when the caller already has them."""
    from scipy import ndimage as ndi

    bad = exp.bad_mask() if bad is None else bad
    sub = exp.background_subtracted() if sub is None else sub
    img = np.where(bad, np.nan, sub).astype(np.float64)
    fin = np.isfinite(img)
    if fin.sum() < 9:
        return np.zeros(img.shape, bool)
    a = img - np.nanmedian(img)
    sig = np.sqrt(np.where(bad | ~(exp.variance > 0), np.nan, exp.variance))
    ring = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    nb_max = ndi.maximum_filter(np.where(fin, a, -np.inf), footprint=ring, mode="constant", cval=-np.inf)
    nb_count = ndi.convolve(fin.astype(int), ring, mode="constant", cval=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return fin & (a > SPIKE_SNR * sig) & (nb_count >= 6) & (nb_max < SPIKE_NEIGHBOUR_RATIO * a)


def reproject_exposure(exp: Exposure, wcs: WCS, shape):
    """Перенос одной экспозиции на общую сетку.

    Возвращает (image, variance, bad, contam): яркость минус фон, дисперсию,
    маску «нельзя использовать» и маску «загрязнение» (блики, засветки — для слоя флагов).
    Плохие родные пиксели (флаги HARD_BAD и одиночные выбросы) заменяются NaN ДО интерполяции:
    билинейная интерполяция не может разнести их значение на соседние пиксели сетки.
    Маски переносятся ближайшим соседом: биты нельзя интерполировать.
    """
    from reproject import reproject_interp

    sub, bad = exp.background_subtracted(), exp.bad_mask()
    bad = bad | lone_spike_mask(exp, sub, bad)
    img = np.where(bad, np.nan, sub).astype(np.float64)
    var = np.where(bad, np.nan, exp.variance).astype(np.float64)
    contam = exp.flag_mask(("CROSSTALK", "GHOST", "GHOST_FPA", "GHOST_EXT", "STREAK",
                            "BLOOM", "SNOWBALL", "HALO", "OVERFLOW", "NONLINEAR"))

    r_img, foot = reproject_interp((img, exp.wcs), wcs, shape_out=shape, order="bilinear")
    r_var, _ = reproject_interp((var, exp.wcs), wcs, shape_out=shape, order="bilinear")
    r_bad, _ = reproject_interp((bad.astype(float), exp.wcs), wcs, shape_out=shape,
                                order="nearest-neighbor")
    r_con, _ = reproject_interp((contam.astype(float), exp.wcs), wcs, shape_out=shape,
                                order="nearest-neighbor")
    outside = foot <= 0
    r_bad = (np.nan_to_num(r_bad, nan=1.0) > 0.5) | outside | ~np.isfinite(r_img) | ~np.isfinite(r_var)
    r_img[r_bad] = np.nan
    r_var[r_bad] = np.nan
    return r_img, r_var, r_bad, np.nan_to_num(r_con) > 0.5


# ---------------------------------------------------------------- кадры
@dataclass
class Frame:
    image: np.ndarray            # сумма эпохи на общей сетке, MJy/sr, фон вычтен; NaN = нет данных
    variance: np.ndarray
    contam: np.ndarray           # bool: загрязнённые/насыщенные пиксели (слой флагов)
    mjd: float                   # средняя дата эпохи
    mjd_min: float
    mjd_max: float
    detector: int
    wave_um: float               # эффективная λ в центре поля
    n_exp: int
    obs_ids: list[str] = field(default_factory=list)
    bg_offsets: list[float | None] = field(default_factory=list)  # остаточный фон каждой экспозиции, MJy/sr
    time_sources: list[dict] = field(default_factory=list)
    wavelength_samples: list[dict] = field(default_factory=list)

    @property
    def span_days(self) -> float:
        return self.mjd_max - self.mjd_min

    @property
    def date(self) -> str:
        return Time(self.mjd, format="mjd").to_value("iso", subfmt="date")


MIN_BG_PIX = 25  # fewer valid pixels: the residual background of an exposure is not estimated
BACKGROUND_METHOD = {
    "name": "zodi_model_plus_clipped_median",
    "description": "Из каждой экспозиции вычитается модель зодиакального света (HDU ZODI), затем остаточный "
                   "постоянный фон: 3σ-клиппированная медиана по пикселям без флага SOURCE на общей сетке.",
    "limitations": [
        "вычитается только постоянная: градиент фона на больших полях остаётся",
        "в плотных полях слабые незафлаганные источники слегка завышают фон",
    ],
}


def _residual_background(img: np.ndarray, good: np.ndarray, source: np.ndarray) -> float | None:
    """Остаточный фон экспозиции. Пиксели с флагом SOURCE исключаются (это сам объект);
    если без них пикселей слишком мало, берутся все хорошие; если и их мало — None (не вычитаем)."""
    for sel in (good & ~source, good):
        if sel.sum() >= MIN_BG_PIX:
            _, med, _ = sigma_clipped_stats(img[sel], sigma=3)
            return float(med)
    return None


def coadd(exps: list[Exposure], wcs: WCS, shape, center: SkyCoord) -> Frame | None:
    """Взвешенное (1/σ²) среднее экспозиций одной эпохи на общей сетке."""
    from reproject import reproject_interp

    num = np.zeros(shape)
    den = np.zeros(shape)
    contam = np.zeros(shape, bool)
    waves, used, offsets, samples = [], [], [], []
    for e in exps:
        img, var, bad, con = reproject_exposure(e, wcs, shape)
        if bad.all():
            continue
        src, _ = reproject_interp((e.flag_mask(("SOURCE",)).astype(float), e.wcs), wcs, shape_out=shape,
                                  order="nearest-neighbor")
        # остаточный фон (модель зодиакального света неидеальна) — убираем медианой по кадру
        med = _residual_background(img, ~bad, np.nan_to_num(src) > 0.5)
        if med is not None:
            img = img - med
        offsets.append(None if med is None else round(med, 6))
        w = np.where(bad, 0.0, 1.0 / np.where(bad, 1.0, var))
        num += np.where(bad, 0.0, img) * w
        den += w
        contam |= con  # включая насыщенные пиксели, которые выброшены из суммы
        if e.contains(center):
            wave, width = e.wavelength_at_sky(center)
            if np.isfinite(wave):
                waves.append(wave)
            samples.append({"obs_id": e.obs_id, "wave_um": float(wave) if np.isfinite(wave) else None,
                            "bw_um": float(width) if np.isfinite(width) else None,
                            "source": "spectral_wcs" if e.wave_wcs is not None else "approximate_detector_position"})
        used.append(e)
    if not used:
        return None
    with np.errstate(invalid="ignore", divide="ignore"):
        image = np.where(den > 0, num / den, np.nan)
        variance = np.where(den > 0, 1.0 / den, np.nan)
    mjds = [e.mjd for e in used]
    return Frame(image=image, variance=variance, contam=contam,
                 mjd=float(np.mean(mjds)), mjd_min=float(min(mjds)), mjd_max=float(max(mjds)),
                 detector=used[0].detector,
                 wave_um=float(np.mean(waves)) if waves else float("nan"),
                 n_exp=len(used), obs_ids=[e.obs_id for e in used], bg_offsets=offsets,
                 time_sources=[{"obs_id": e.obs_id, "source": e.time_source,
                                "assumptions": e.time_assumptions} for e in used],
                 wavelength_samples=samples)


def spectral_sampling(frames: list[Frame]) -> dict:
    """Which wavelengths the target was seen at, per exposure (LVF: each exposure has its own λ at a point).

    different_channels=True when the spread of λ is larger than one spectral channel (the median width of
    the samples, about λ/R). Then brightness differences between epochs mix variability with the spectrum."""
    samples = [s for f in frames for s in f.wavelength_samples if s.get("wave_um") is not None]
    if not samples:
        return {"wave_min_um": None, "wave_max_um": None, "channel_width_um": None, "different_channels": None,
                "n_samples": 0}
    w = np.array([s["wave_um"] for s in samples], float)
    bw = np.array([s["bw_um"] for s in samples if s.get("bw_um")], float)
    width = float(np.median(bw)) if bw.size else None
    spread = float(w.max() - w.min())
    limit = width if width else 0.05 * float(np.mean(w))
    return {"wave_min_um": round(float(w.min()), 4), "wave_max_um": round(float(w.max()), 4),
            "channel_width_um": None if width is None else round(width, 4),
            "different_channels": bool(spread > limit), "n_samples": len(samples)}


MIN_DIFF_PIX = 50  # fewer valid pixels: the empirical noise of a difference frame is not estimated
DIFF_METHOD = {
    "name": "leave_one_out_median",
    "description": "Кадр минус медиана ОСТАЛЬНЫХ эпох (без самого кадра), делённая на шум; шум затем "
                   "нормирован робастно (MAD) по самому разностному кадру. Шкала показывает относительные остатки.",
    "units": "robust_residual_scale",
    "limitations": [
        "это не откалиброванная статистическая значимость: корреляции интерполяции и оценка фона не распространены полностью",
        "λ в точке может различаться между эпохами (LVF): у объектов с крутым спектром остаток бывает не движением",
        "переменная звезда даёт остаток того же знака на месте, а не пару «+/−»",
        "медленное движение (меньше ~1 пикселя между эпохами) даёт слабую пару, слившуюся в один диполь",
    ],
}


def _robust_std(d: np.ndarray) -> float | None:
    v = d[np.isfinite(d)]
    if v.size < MIN_DIFF_PIX:
        return None
    s = float(1.4826 * np.median(np.abs(v - np.median(v))))
    return s if s > 0 else None


def difference_details(frames: list[Frame], renorm: bool = True) -> dict:
    """Разностные кадры в единицах σ и сведения о нормировке.

    Шаблон для кадра i — медиана ОСТАЛЬНЫХ эпох: если включать в него сам кадр, при двух эпохах
    шаблон равен среднему и остаток вдвое слабее, а шум кадра и шаблона коррелирован.
    Дисперсия медианы m значений: σ²/m для m ≤ 2 (это среднее), ≈ (π/2)·σ²/m для m ≥ 3.
    Распространённая дисперсия после билинейного перепроецирования завышена (интерполяция сглаживает
    шум), поэтому при renorm результат делится на робастное σ самого разностного кадра (noise_scale).
    Движущийся объект даёт пару: «+» там, где он в этом кадре, «−» там, где он в остальных.
    -> {diffs, template (медиана всех эпох, фон для «следа»), noise_scale: [float|None], method}
    """
    import warnings

    cube = np.array([f.image for f in frames])
    vcube = np.array([f.variance for f in frames])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        template = np.nanmedian(cube, axis=0)
        diffs, scales = [], []
        for i, f in enumerate(frames):
            if len(frames) < 2:
                diffs.append(np.where(np.isfinite(f.image), 0.0, np.nan))
                scales.append(None)
                continue
            others, ovar = np.delete(cube, i, 0), np.delete(vcube, i, 0)
            ref = np.nanmedian(others, axis=0)
            m = np.isfinite(others).sum(axis=0)
            factor = np.where(m <= 2, 1.0, np.pi / 2) / np.maximum(m, 1)
            var_ref = np.nanmean(ovar, axis=0) * factor
            with np.errstate(invalid="ignore", divide="ignore"):
                d = (f.image - ref) / np.sqrt(f.variance + var_ref)
            d[~np.isfinite(d)] = np.nan
            sc = _robust_std(d) if renorm else None
            if sc is not None:
                d = d / sc
            diffs.append(d)
            scales.append(None if sc is None else round(sc, 4))
    return {"diffs": diffs, "template": template, "noise_scale": scales, "method": DIFF_METHOD}


def difference(frames: list[Frame], renorm: bool = True) -> tuple[list[np.ndarray], np.ndarray]:
    """Разностные кадры в единицах σ и шаблон (см. difference_details)."""
    r = difference_details(frames, renorm)
    return r["diffs"], r["template"]


TRAIL_OLD = (0.25, 0.55, 1.0)   # синий — самая ранняя эпоха
TRAIL_NEW = (1.0, 0.85, 0.15)   # жёлтый — самая поздняя


def trail(diff_sigma: list[np.ndarray], clip: float = 8.0, floor: float = 2.0) -> np.ndarray:
    """«След» движения на одной картинке (RGB 0..1).

    Положительные остатки каждой эпохи окрашиваются цветом по времени (старые — синие,
    средние — зелёные, новые — жёлтые) с прозрачностью α от 0.55 до 1.0 и СКЛАДЫВАЮТСЯ.
    Обычное наложение слоёв не годится: слой с α=1 полностью закрыл бы все предыдущие.
    Все цвета палитры яркие: у viridis начало почти чёрное, и старая эпоха терялась на тёмном фоне.
    """
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("time", [TRAIL_OLD, (0.2, 0.9, 0.55), TRAIL_NEW])
    n = len(diff_sigma)
    rgb = np.zeros((*diff_sigma[0].shape, 3))
    for i, d in enumerate(diff_sigma):
        t = i / max(n - 1, 1)
        alpha = 0.55 + 0.45 * t
        # остатки ниже floor σ — шум: без порога каждая эпоха даёт цветную «рябь» по всему кадру
        pos = np.clip(np.nan_to_num(d) - floor, 0, clip - floor) / (clip - floor)
        rgb += alpha * pos[..., None] * np.array(cmap(t)[:3])
    return np.clip(rgb, 0, 1)
