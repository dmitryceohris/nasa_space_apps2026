"""Measurements: centroid, the target's track over epochs, proper-motion fit, spectrum and light curve."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time

from Exposure import Exposure
from processing.centroid import METHOD as CENTROID_METHOD
from processing.centroid import measure as _measure
from processing.photometry import PHOTOMETRY_METHOD
from processing.photometry import light_curve as _light_curve
from processing.photometry import spectrum as _spectrum

from .data import target as _target
from .stack import Stack

ASSOC_PX = 1.0  # Same conservative display gate as /api/view, not a calibrated association radius.
# Agreement with the prediction does not establish target identity or registration accuracy.
ASSIST_PX = 3.0  # with assist=True: how far from the prediction a candidate may be picked or a window started
ASSISTED = ("nearest_to_prediction", "windowed_at_prediction")
WINDOW_R = 2.5  # radius (px) of the windowed centroid used at the predicted position
MIN_FORCED_SNR = 5.0  # a windowed measurement below this S/N is not accepted


def centroid(img, signed: bool = False, saturated=None) -> dict:
    """Centre of the object in any 2-D array (numpy pixels, y up).
    -> {status 'ok'|'ambiguous'|'none', reason, xy, err_px, quality 'good'|'fair'|'poor', flags, candidates, ...}
    signed=True for difference frames: '+' and '-' residuals separately, the position comes from '+'."""
    return _measure(np.asarray(img, float), signed, saturated=saturated)


def windowed_centroid(img, x0: float, y0: float, r: float = WINDOW_R, iters: int = 20) -> dict | None:
    """Iterative flux-weighted centroid inside a circle of radius r, starting at (x0, y0) (like SExtractor XWIN).
    Background and noise come from the whole frame (median / MAD). Used where the target is known to be:
    it does not merge the target with a neighbour whose wings touch it. -> {xy, err_px, snr} or None."""
    from processing.centroid import _robust

    a = np.asarray(img, float)
    fin = np.isfinite(a)
    _, sig = _robust(np.where(fin, a, np.nan))
    if not sig > 0:
        return None
    yy, xx = np.indices(a.shape)

    def local_bg(x, y):
        """Plane a + b*dx + c*dy fitted (3-sigma clipped) in the annulus r+1..r+4 around (x, y): removes the
        sloping wing of a neighbour, which would otherwise pull the centre towards it."""
        rr = np.hypot(xx - x, yy - y)
        ann = fin & (rr >= r + 1) & (rr <= r + 4)
        if ann.sum() < 8:
            return np.full(a.shape, np.nanmedian(a[fin]))
        A = np.column_stack([np.ones(ann.sum()), xx[ann] - x, yy[ann] - y])
        v, keep = a[ann], np.ones(ann.sum(), bool)
        for _ in range(3):
            coef, *_ = np.linalg.lstsq(A[keep], v[keep], rcond=None)
            res = v - A @ coef
            keep = np.abs(res) < 3 * max(sig, 1.4826 * np.median(np.abs(res[keep])))
            if keep.sum() < 8:
                break
        return coef[0] + coef[1] * (xx - x) + coef[2] * (yy - y)

    x, y = float(x0), float(y0)
    for _ in range(iters):
        bg = local_bg(x, y)
        m = fin & (np.hypot(xx - x, yy - y) <= r)
        w = np.where(m, np.clip(a - bg, 0, None), 0.0)
        t = w.sum()
        if t <= 0:
            return None
        nx, ny = float((w * xx).sum() / t), float((w * yy).sum() / t)
        moved = np.hypot(nx - x, ny - y)
        x, y = nx, ny
        if moved < 1e-3:
            break
    bg = local_bg(x, y)
    m = fin & (np.hypot(xx - x, yy - y) <= r)
    w = np.where(m, np.clip(a - bg, 0, None), 0.0)
    t = w.sum()
    if t <= 0 or m.sum() == 0:
        return None
    snr = float(np.sum((a - bg)[m]) / (sig * np.sqrt(m.sum())))
    ex = sig * np.sqrt(((xx[m] - x) ** 2).sum()) / t
    ey = sig * np.sqrt(((yy[m] - y) ** 2).sum()) / t
    return {"xy": [x, y], "err_px": [float(ex), float(ey)], "snr": snr}


# ------------------------------------------------------------------ annual parallax
def parallax_offset(ra_deg, dec_deg, mjd, parallax_mas) -> np.ndarray:
    """Shift (arcsec, [east, north]) of the geocentric position relative to the barycentric one.
    Standard parallax factors with the Earth's barycentric position (astropy built-in ephemeris).
    SPHEREx is in low Earth orbit: its own offset from the geocentre is ~1e-5 of this, neglected."""
    if not parallax_mas:
        return np.zeros((np.size(mjd), 2)).squeeze()
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    e = get_body_barycentric("earth", t).xyz.to_value(u.au)  # (3, N), ICRS axes
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    p = parallax_mas / 1000
    east = p * (e[0] * np.sin(a) - e[1] * np.cos(a))
    north = p * (e[0] * np.cos(a) * np.sin(d) + e[1] * np.sin(a) * np.sin(d) - e[2] * np.cos(d))
    out = np.stack([east, north], axis=-1)
    return out[0] if np.ndim(mjd) == 0 else out


def predict(where, mjd) -> SkyCoord:
    """Catalogue position on `mjd` seen from the Earth: space motion + annual parallax."""
    t = _target(where)
    c = t.at(float(mjd))
    if t.parallax > 0:
        de, dn = parallax_offset(c.ra.deg, c.dec.deg, float(mjd), t.parallax)
        c = c.spherical_offsets_by(de * u.arcsec, dn * u.arcsec)
    return c


# ------------------------------------------------------------------ track
TRACK_COLUMNS = ("mjd", "x", "y", "err_x", "err_y", "ra", "dec", "err_arcsec", "pred_x", "pred_y", "offset_px")


@dataclass
class Track:
    """Position of the target in every epoch. rows: list of dicts; to_numpy() -> (N, len(columns)) floats,
    NaN where nothing was measured. Columns: mjd, x, y, err_x, err_y (px), ra, dec (deg), err_arcsec,
    pred_x, pred_y (px), offset_px (measured minus predicted)."""
    rows: list[dict]
    target: object
    pixel_scale: float
    diff: bool = False

    columns = TRACK_COLUMNS

    def to_numpy(self, columns=TRACK_COLUMNS) -> np.ndarray:
        return np.array([[np.nan if r.get(c) is None else r[c] for c in columns] for r in self.rows], float) \
            .reshape(len(self.rows), len(columns))

    def __array__(self, dtype=None, copy=None):
        return np.asarray(self.to_numpy(), dtype=dtype)

    def good(self, include_poor: bool = False, include_assisted: bool = False) -> list[dict]:
        """Rows usable for a motion fit: measured, near the prediction (if any), quality good/fair.
        Rows found with the help of the prediction (assist=True) only when include_assisted=True."""
        ok = []
        for r in self.rows:
            if r["xy"] is None or "far_from_prediction" in r["flags"]:
                continue
            if r["quality"] == "poor" and not include_poor:
                continue
            if r["association"] in ASSISTED and not include_assisted:
                continue
            ok.append(r)
        return ok

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, i):
        return self.rows[i]

    def __repr__(self) -> str:
        n = sum(r["xy"] is not None for r in self.rows)
        return f"Track({self.target.name!r}, {n}/{len(self.rows)} epochs measured)"


def track(st: Stack, diff: bool = False, assoc_px: float = ASSOC_PX, assist: bool = False) -> Track:
    """Measures the target in every epoch of the stack.

    Default (same rule as /api/view): the position is taken only when the frame has one unambiguous object
    (association='unique'); farther than assoc_px (1 px) from the prediction it is flagged
    'far_from_prediction' and not used by fit_motion. The prediction never chooses the object.

    assist=True (exploration, NOT an independent measurement): the prediction may help, within ASSIST_PX:
      * ambiguous frame -> the detected candidate nearest to the prediction ('nearest_to_prediction');
      * nothing usable near the prediction (no detection, or the detection is a blend centred elsewhere) ->
        windowed centroid started at the prediction, accepted at S/N >= 5
        ('windowed_at_prediction', flag 'blended' or 'forced', quality at best 'fair').
    Such rows are used by fit_motion only with include_assisted=True: the prediction chose them, so a motion
    fitted from them is not independent of the catalogue.
    diff=True measures the '+' residual of the difference frames instead of the coadds."""
    from .stack import diff as _diff

    imgs = _diff(st) if diff else st.cube
    has_model = bool(st.target.pmra or st.target.pmdec or st.target.parallax)
    rows = []
    for i, f in enumerate(st.frames):
        m = centroid(imgs[i], signed=diff, saturated=f.contam)
        pred = None
        if has_model:
            px, py = st.pixel(predict(st.target, f.mjd))
            pred = [float(px), float(py)]
        xy, err, quality, flags, assoc = m["xy"], m["err_px"], m["quality"], list(m["flags"]), None
        if xy is not None:
            assoc = "unique"
        elif assist and pred is not None and m["status"] == "ambiguous":
            cands = [c for c in m["candidates"] if c["sign"] in (None, 1)]
            if cands:
                c = min(cands, key=lambda c: np.hypot(c["xy"][0] - pred[0], c["xy"][1] - pred[1]))
                if np.hypot(c["xy"][0] - pred[0], c["xy"][1] - pred[1]) <= ASSIST_PX:
                    xy, err, quality, flags = c["xy"], c["err_px"], c["quality"], list(c["flags"])
                    assoc = "nearest_to_prediction"
        limit = ASSIST_PX if assoc == "nearest_to_prediction" else assoc_px
        far = xy is not None and pred is not None and np.hypot(xy[0] - pred[0], xy[1] - pred[1]) > limit
        if assist and pred is not None and (xy is None or far):
            w = windowed_centroid(imgs[i], *pred)  # weights are clipped at 0: in diff only the "+" counts
            if (w is not None and w["snr"] >= MIN_FORCED_SNR
                    and np.hypot(w["xy"][0] - pred[0], w["xy"][1] - pred[1]) <= ASSIST_PX):
                flags = ["blended"] if far else ["forced"]
                xy, err, quality, assoc = w["xy"], w["err_px"], "fair", "windowed_at_prediction"
        row = {"index": i, "date": f.date, "mjd": f.mjd, "wave_um": f.wave_um, "n_exp": f.n_exp,
               "status": m["status"], "reason": m["reason"], "association": assoc,
               "xy": xy, "x": None, "y": None, "err_x": None, "err_y": None, "ra": None, "dec": None,
               "err_arcsec": None, "quality": quality if xy is not None else None, "flags": flags if xy else [],
               "pred_xy": pred, "pred_x": pred[0] if pred else None, "pred_y": pred[1] if pred else None,
               "offset_px": None}
        if xy is not None:
            ra, dec = st.radec(*xy)
            row.update(x=xy[0], y=xy[1], err_x=err[0], err_y=err[1], ra=float(ra), dec=float(dec),
                       err_arcsec=float(np.hypot(*err) / np.sqrt(2) * st.pixel_scale))
            if pred is not None:
                off = float(np.hypot(xy[0] - pred[0], xy[1] - pred[1]))
                row["offset_px"] = off
                if off > (ASSIST_PX if assoc in ASSISTED else assoc_px):
                    row["flags"] = row["flags"] + ["far_from_prediction"]
        rows.append(row)
    return Track(rows, st.target, st.pixel_scale, diff)


# ------------------------------------------------------------------ proper motion
@dataclass
class Motion:
    """Measured proper motion. pm in arcsec/yr: [mu_ra*cos(dec), mu_dec] with 1-sigma errors.
    to_numpy() -> array([pm_ra, pm_dec]); to_dict() -> everything."""
    pm: np.ndarray
    pm_err: np.ndarray
    total: float
    total_err: float
    position_angle_deg: float
    n_epochs: int
    span_years: float
    chi2_red: float | None
    errors_scaled: bool
    parallax_mode: str
    parallax_mas: float | None
    reference_pm: np.ndarray | None
    difference_sigma: np.ndarray | None
    residuals_arcsec: np.ndarray
    used_indices: list[int]
    warnings: list[str] = field(default_factory=list)
    method: str = ("взвешенный МНК по положениям эпох: позиция + собственное движение "
                   "(+ годичный параллакс, если задан или подгоняется)")

    def to_numpy(self) -> np.ndarray:
        return np.asarray(self.pm, float)

    def __array__(self, dtype=None, copy=None):
        return np.asarray(self.pm, dtype=dtype)

    def to_dict(self) -> dict:
        out = {}
        for k, v in self.__dict__.items():
            out[k] = v.tolist() if isinstance(v, np.ndarray) else v
        return out

    def __repr__(self) -> str:
        s = (f"Motion(pm = [{self.pm[0]:+.2f} ± {self.pm_err[0]:.2f}, {self.pm[1]:+.2f} ± {self.pm_err[1]:.2f}]"
             f" arcsec/yr, {self.n_epochs} epochs over {self.span_years:.2f} yr")
        if self.reference_pm is not None:
            s += f"; catalogue [{self.reference_pm[0]:+.2f}, {self.reference_pm[1]:+.2f}]"
        return s + ")"


def fit_motion(tr: Track, parallax: str = "auto", include_poor: bool = False,
               include_assisted: bool = False) -> Motion:
    """Proper motion from the measured positions of a Track.

    Model per axis: offset(t) = p0 + mu*(t - t_mean) + plx * P(t), P = parallax factors.
    parallax: 'auto' = fixed to the catalogue value when the target has one, else ignored;
              'catalog' / 'none' / 'fit' (needs >= 4 epochs).
    Errors are the centroid noise errors (in arcsec); when chi2_red > 1 they are scaled by sqrt(chi2_red),
    because the noise error does not include PSF or registration systematics."""
    rows = tr.good(include_poor, include_assisted)
    warn = []
    if len(rows) < 2:
        raise ValueError(f"нужно минимум 2 измеренные эпохи, есть {len(rows)}")
    mjd = np.array([r["mjd"] for r in rows])
    ra = np.array([r["ra"] for r in rows])
    dec = np.array([r["dec"] for r in rows])
    sig = np.array([max(r["err_arcsec"], 1e-3) for r in rows])
    # mean RA through the 0/360 wrap: 359.99 and 0.01 average to 0.0, not 180
    ra0 = (ra[0] + np.mean((ra - ra[0] + 180) % 360 - 180)) % 360
    dec0 = np.mean(dec)
    east = (ra - ra0 + 180) % 360 - 180
    east = east * np.cos(np.radians(dec0)) * 3600
    north = (dec - dec0) * 3600
    t_yr = (mjd - mjd.mean()) / 365.25
    span = float(np.ptp(mjd) / 365.25)

    plx_known = float(getattr(tr.target, "parallax", 0.0) or 0.0)
    mode = parallax
    if mode == "auto":
        mode = "catalog" if plx_known > 0 else "none"
    if mode == "catalog" and plx_known <= 0:
        raise ValueError("у цели нет параллакса в каталоге: используйте parallax='none' или 'fit'")
    if mode == "fit" and len(rows) < 4:
        raise ValueError("подгонка параллакса требует минимум 4 эпох")
    if mode not in ("catalog", "none", "fit"):
        raise ValueError("parallax: 'auto', 'catalog', 'none' или 'fit'")

    pf = parallax_offset(ra0, dec0, mjd, 1000.0)  # (N, 2) shift for a 1-arcsec parallax
    pf = pf.reshape(len(rows), 2)
    if mode == "catalog":
        east = east - pf[:, 0] * plx_known / 1000
        north = north - pf[:, 1] * plx_known / 1000

    # joint weighted least squares on both axes: unknowns [e0, n0, mu_e, mu_n, (plx)]
    n = len(rows)
    ncol = 5 if mode == "fit" else 4
    A = np.zeros((2 * n, ncol))
    A[:n, 0], A[n:, 1] = 1, 1
    A[:n, 2], A[n:, 3] = t_yr, t_yr
    if mode == "fit":
        A[:n, 4], A[n:, 4] = pf[:, 0], pf[:, 1]
    b = np.concatenate([east, north])
    w = np.concatenate([1 / sig, 1 / sig])
    Aw, bw = A * w[:, None], b * w
    sol, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    cov = np.linalg.pinv(Aw.T @ Aw)
    res = b - A @ sol
    dof = 2 * n - ncol
    chi2_red = float(np.sum((res * w) ** 2) / dof) if dof > 0 else None
    scaled = bool(chi2_red is not None and chi2_red > 1)
    if scaled:
        cov = cov * chi2_red
    if dof <= 0:
        warn.append("степеней свободы нет: ошибки только шумовые и, вероятно, занижены")
    if span < 0.25:
        warn.append(f"интервал всего {span:.2f} года: собственное движение определено грубо")
    pm = sol[2:4]
    pm_err = np.sqrt(np.diag(cov)[2:4])
    total = float(np.hypot(*pm))
    total_err = float(np.sqrt((pm[0] * pm_err[0]) ** 2 + (pm[1] * pm_err[1]) ** 2) / total) if total > 0 else 0.0
    pa = float(np.degrees(np.arctan2(pm[0], pm[1])) % 360)  # from north through east
    ref = diff = None
    if getattr(tr.target, "pmra", 0) or getattr(tr.target, "pmdec", 0):
        ref = np.array([tr.target.pmra, tr.target.pmdec]) / 1000
        diff = (pm - ref) / np.where(pm_err > 0, pm_err, np.nan)
    return Motion(pm=pm, pm_err=pm_err, total=total, total_err=total_err, position_angle_deg=pa, n_epochs=n,
                  span_years=span, chi2_red=chi2_red, errors_scaled=scaled, parallax_mode=mode,
                  parallax_mas=float(sol[4] * 1000) if mode == "fit" else (plx_known if mode == "catalog" else None),
                  reference_pm=ref, difference_sigma=diff,
                  residuals_arcsec=np.stack([res[:n], res[n:]], axis=1), used_indices=[r["index"] for r in rows],
                  warnings=warn)


# ------------------------------------------------------------------ spectrum and light curve
class _Points:
    """List of photometry points with numpy views. to_numpy() -> (N, 3) float."""

    def __init__(self, points: list[dict]):
        self.points = points

    def _arr(self, key, only_ok):
        pts = [p for p in self.points if p.get("ok") or not only_ok]
        return np.array([np.nan if p.get(key) is None else p[key] for p in pts], float)

    @property
    def ok(self) -> np.ndarray:
        return np.array([bool(p.get("ok")) for p in self.points])

    def __len__(self):
        return len(self.points)

    def __iter__(self):
        return iter(self.points)

    def __getitem__(self, i):
        return self.points[i]

    def __array__(self, dtype=None, copy=None):
        return np.asarray(self.to_numpy(), dtype=dtype)


class Spectrum(_Points):
    """Spectrum points: wave_um, flux_ujy, err_ujy, ok (False = bad/saturated pixels in the aperture), det, mjd."""
    method = PHOTOMETRY_METHOD

    def to_numpy(self, only_ok: bool = True) -> np.ndarray:
        """(N, 3): wavelength um, flux uJy, error uJy; only reliable points unless only_ok=False."""
        return np.stack([self._arr("wave_um", only_ok), self._arr("flux_ujy", only_ok),
                         self._arr("err_ujy", only_ok)], axis=1).reshape(-1, 3)

    @property
    def wave(self):
        return self._arr("wave_um", False)

    @property
    def flux(self):
        return self._arr("flux_ujy", False)

    @property
    def err(self):
        return self._arr("err_ujy", False)

    def __repr__(self):
        return f"Spectrum({len(self)} points, {int(self.ok.sum())} reliable)"


class LightCurve(_Points):
    """Flux over time in a wavelength window."""

    def __init__(self, points, wave_min, wave_max):
        super().__init__(points)
        self.window_um = (wave_min, wave_max)

    def to_numpy(self, only_ok: bool = True) -> np.ndarray:
        """(N, 3): mjd, flux uJy, error uJy."""
        return np.stack([self._arr("mjd", only_ok), self._arr("flux_ujy", only_ok),
                         self._arr("err_ujy", only_ok)], axis=1).reshape(-1, 3)

    @property
    def mjd(self):
        return self._arr("mjd", False)

    @property
    def flux(self):
        return self._arr("flux_ujy", False)

    def __repr__(self):
        return f"LightCurve({len(self)} points in {self.window_um[0]}-{self.window_um[1]} um)"


def spectrum(exps: list[Exposure], where) -> Spectrum:
    """Forced aperture photometry at the target (moved to each exposure's date) on the native pixels of every
    exposure; each exposure sees the target at its own wavelength (LVF), together they form the spectrum."""
    t = _target(where)
    return Spectrum(_spectrum(exps, lambda mjd: predict(t, mjd)))


def lightcurve(spec: Spectrum, wave_min: float, wave_max: float) -> LightCurve:
    """Points of the spectrum inside [wave_min, wave_max] um, sorted by time. A narrow window keeps the
    spectral shape from being mistaken for variability."""
    return LightCurve(_light_curve(spec.points, wave_min, wave_max), wave_min, wave_max)
