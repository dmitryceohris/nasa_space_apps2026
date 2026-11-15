"""Epoch stacks on one sky grid, difference frames, trail, PNG/GIF/FITS output."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from Exposure import Exposure
from processing.epochs import GAP_DAYS, MAX_SPAN_DAYS
from processing.frames import Frame, difference_details
from processing.frames import trail as _trail
from processing.pipeline import frames_from_exposures, render_view
from processing.targets import Target

from .data import target as _target


@dataclass
class Stack:
    """Coadded epochs on one north-up TAN grid (native SPHEREx pixel, 6.15").

    Arrays are numpy with shape (N_epochs, H, W); pixel [x, y] is 0-based, y UP (row 0 = south), east left.
    np.asarray(stack) and stack.to_numpy() give the image cube."""
    frames: list[Frame]
    wcs: object
    target: Target
    params: dict = field(default_factory=dict)
    exposures: list[Exposure] = field(default_factory=list,repr=False)

    # ---- arrays
    @cached_property
    def cube(self) -> np.ndarray:
        """(N, H, W) MJy/sr, background removed, NaN = no data."""
        return np.stack([f.image for f in self.frames]) if self.frames else np.empty((0, *self.shape))

    @cached_property
    def variance(self) -> np.ndarray:
        return np.stack([f.variance for f in self.frames]) if self.frames else np.empty((0, *self.shape))

    @cached_property
    def flags(self) -> np.ndarray:
        """(N, H, W) bool: contaminated or saturated pixels (ghosts, streaks, overflow, ...)."""
        return np.stack([f.contam for f in self.frames]) if self.frames else np.empty((0, *self.shape), bool)

    @property
    def mjd(self) -> np.ndarray:
        return np.array([f.mjd for f in self.frames])

    @property
    def dates(self) -> list[str]:
        return [f.date for f in self.frames]

    @property
    def wave(self) -> np.ndarray:
        """Wavelength (um) at the target in each epoch; NaN when unknown."""
        return np.array([f.wave_um for f in self.frames])

    @property
    def n_exp(self) -> np.ndarray:
        return np.array([f.n_exp for f in self.frames])

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(self.wcs.pixel_shape[::-1]) if self.wcs.pixel_shape else self.frames[0].image.shape

    @property
    def pixel_scale(self) -> float:
        """Arcsec per pixel."""
        return float(abs(self.wcs.wcs.cdelt[1]) * 3600)

    def to_numpy(self) -> np.ndarray:
        return self.cube

    def __array__(self, dtype=None, copy=None):
        return np.asarray(self.cube, dtype=dtype)

    def __len__(self) -> int:
        return len(self.frames)

    def __repr__(self) -> str:
        span = f"{self.dates[0]}..{self.dates[-1]}" if self.frames else "empty"
        return f"Stack({self.target.name!r}, {len(self)} epochs {span}, {self.shape[1]}x{self.shape[0]} px)"

    # ---- coordinates
    def radec(self, x, y) -> tuple:
        """numpy pixel -> (ra_deg, dec_deg)."""
        c = self.wcs.pixel_to_world(x, y)
        return c.ra.deg, c.dec.deg

    def pixel(self, ra, dec=None) -> tuple:
        """(ra_deg, dec_deg) or SkyCoord -> numpy pixel (x, y)."""
        c = ra if isinstance(ra, SkyCoord) else SkyCoord(ra, dec, unit="deg")
        x, y = self.wcs.world_to_pixel(c)
        return x, y

    # ---- differences (computed once)
    @cached_property
    def _diff(self) -> dict | None:
        return difference_details(self.frames) if len(self.frames) > 1 else None


def stack(exps: list[Exposure], where=None, fov: float = 3.0, gap: float = GAP_DAYS,
          max_span: float = MAX_SPAN_DAYS) -> Stack:
    """Groups exposures into epochs and coadds each epoch (1/variance weights, zodi and residual background
    removed, bad pixels masked) on one grid centred on the target at the mean date of the exposures.
    where: name / (ra, dec) / SkyCoord / Target; default: centre of the first exposure."""
    if not exps:
        raise ValueError("нет снимков")
    if where is None:
        e = exps[0]
        where = e.pixel_to_sky((e.shape[1] - 1) / 2, (e.shape[0] - 1) / 2)
    t = _target(where)
    mjd = float(np.mean([e.mjd for e in exps]))
    frames, wcs, shape = frames_from_exposures(exps, t.at(mjd), fov, gap, max_span)
    wcs.pixel_shape = (shape[1], shape[0])
    return Stack(frames, wcs, t, {"fov_arcmin": fov, "gap_days": gap, "max_span_days": max_span,
                                  "n_exposures": len(exps), "center_mjd": mjd},exposures=list(exps))


def diff(st: Stack, info: bool = False):
    """(N, H, W) difference frames in robust residual units, not calibrated statistical sigma.
    Each epoch minus the median of the OTHER epochs. Positive = brighter now.
    info=True -> (cube, {noise_scale, template, method})."""
    if st._diff is None:
        raise ValueError("для разности нужно минимум 2 эпохи")
    cube = np.stack(st._diff["diffs"])
    if not info:
        return cube
    return cube, {"noise_scale": st._diff["noise_scale"], "template": st._diff["template"],
                  "method": st._diff["method"]}


def trail(st: Stack) -> np.ndarray:
    """(H, W, 3) RGB 0..1: positive residuals of every epoch coloured by time (blue first, yellow last)."""
    if st._diff is None:
        raise ValueError("для следа нужно минимум 2 эпохи")
    return _trail(st._diff["diffs"])


def save(st: Stack, out_dir, url_prefix: str = ".") -> dict:
    """PNG of every epoch (coadd, diff, flags), trail, GIF and result.json, as the web site makes them.
    PNG are north up (y flipped). -> the JSON."""
    out = Path(out_dir)
    res = render_view(st.frames, st.wcs, st.shape, out, url_prefix, st.target)
    (out / "result.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    return res


def to_fits(st: Stack, path, what: str = "coadd", overwrite: bool = False) -> Path:
    """FITS cube with the celestial WCS, for DS9 / Aladin / astropy.
    HDUs: PRIMARY image or robust normalized difference, native coadd VARIANCE, FLAGS,
    EPOCHS table (index, mjd, mjd_min, mjd_max, wave_um, n_exp)."""
    from astropy.io import fits
    from astropy.table import Table

    if what not in ("coadd", "diff"):
        raise ValueError("what: 'coadd' или 'diff'")
    data = st.cube if what == "coadd" else diff(st)
    hdr = st.wcs.to_header()
    hdr["BUNIT"] = "MJy/sr" if what == "coadd" else "1"
    if what == "diff":
        hdr["DIFFUNIT"] = "robust_residual_scale"
        hdr["HISTORY"] = "Not calibrated statistical significance; VARIANCE is the coadd variance."
    hdr["CONTENT"] = (what, "coadd image or diff vs other epochs")
    hdr["OBJECT"] = st.target.name
    hdr["ORIGIN"] = "spxhack (NASA SPHEREx via NASA/IPAC IRSA)"
    hdr["DATE"] = Time.now().isot[:19]
    prim = fits.PrimaryHDU(data.astype(np.float32), header=hdr)
    vh = st.wcs.to_header()
    vh["BUNIT"] = "(MJy/sr)^2"
    vh["VARTYPE"] = ("COADD", "Variance of coadd, not normalized difference")
    var = fits.ImageHDU(st.variance.astype(np.float32), header=vh, name="VARIANCE")
    flg = fits.ImageHDU(st.flags.astype(np.uint8), header=st.wcs.to_header(), name="FLAGS")
    tab = Table({"index": np.arange(len(st)), "mjd": st.mjd,
                 "mjd_min": [f.mjd_min for f in st.frames], "mjd_max": [f.mjd_max for f in st.frames],
                 "wave_um": st.wave, "n_exp": st.n_exp})
    ep = fits.table_to_hdu(tab)
    ep.name = "EPOCHS"
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([prim, var, flg, ep]).writeto(p, overwrite=overwrite)
    return p
