"""spxhack: SPHEREx images over time, in your own Python code.

    import spxhack as spx

    t    = spx.target("barnard")                 # name, (ra, dec) or SkyCoord; t.at(mjd) -> position on a date
    rows = spx.search(t, det=3)                  # which exposures exist (network, no download)
    exps = spx.fetch(t, det=3, fov=3)            # download cutouts (network) ...
    exps = spx.load(["a.fits", "b.fits"])        # ... or open local files / URLs

    st  = spx.stack(exps, t, fov=3)              # epochs -> coadds on one north-up grid
    st.cube                                      # np.ndarray (N, H, W), MJy/sr; also np.asarray(st), st.to_numpy()
    d   = spx.diff(st)                           # (N, H, W) in sigma
    rgb = spx.trail(st)                          # (H, W, 3)

    spx.centroid(img)                            # centre on any 2-D array: xy, err, quality, status
    tr = spx.track(st)                           # the target in every epoch; tr.to_numpy() -> (N, 11)
    pm = spx.fit_motion(tr)                      # measured proper motion vs catalogue; pm.to_numpy()
    sp = spx.spectrum(exps, t)                   # photometry at every wavelength; sp.to_numpy() -> (N, 3)
    lc = spx.lightcurve(sp, 2.0, 2.2)            # flux over time in a window; lc.to_numpy() -> (N, 3)

    spx.show(st)                                 # blinking animation (inline in Jupyter)
    spx.save(st, "out/")                         # PNG, GIF, result.json as on the web site
    spx.to_fits(st, "cube.fits")                 # cube + WCS for DS9 / Aladin

Pixels are numpy ones: [x, y], 0-based, y grows to the north, x to the west (east is left). Every measurement
carries its method, error and quality; see spx.METHODS. Not done: PSF fitting, deblending, neural networks.
"""
from Exposure import Exposure
from processing.centroid import METHOD as _CENTROID
from processing.frames import BACKGROUND_METHOD as _BACKGROUND
from processing.frames import DIFF_METHOD as _DIFF
from processing.photometry import PHOTOMETRY_METHOD as _PHOT
from processing.targets import Target

from .data import epochs, fetch, load, search, target
from .measure import (LightCurve, Motion, Spectrum, Track, centroid, fit_motion, lightcurve, parallax_offset,
                      predict, spectrum, track, windowed_centroid)
from .show import show
from .stack import Stack, diff, save, stack, to_fits, trail

__version__ = "0.1.0"

METHODS = {"background": _BACKGROUND, "difference": _DIFF, "centroid": _CENTROID, "photometry": _PHOT}

__all__ = [
    "target", "search", "fetch", "load", "epochs",
    "stack", "diff", "trail", "save", "to_fits",
    "centroid", "windowed_centroid", "track", "fit_motion", "predict", "parallax_offset", "spectrum", "lightcurve",
    "show",
    "Exposure", "Target", "Stack", "Track", "Motion", "Spectrum", "LightCurve", "METHODS", "__version__",
]
