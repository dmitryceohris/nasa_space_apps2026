"""numpy -> PNG. The vertical flip lives ONLY here: PNG has y down, FITS/numpy has y up (north on top after flip).
Every pixel coordinate sent to the front goes through png_xy(), otherwise markers drift off the objects."""
import numpy as np
from astropy.wcs import WCS
from PIL import Image

PIX_ARCSEC = 6.15  # SPHEREx pixel scale (s_pixel_scale in the IRSA exposure table)


def grid(ra, dec, fov_arcmin):
    """Common TAN grid centred on the target; -> (WCS, (ny, nx))."""
    n = int(np.clip(round(fov_arcmin * 60 / PIX_ARCSEC), 8, 400))
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-PIX_ARCSEC / 3600, PIX_ARCSEC / 3600]
    return w, (n, n)


def wcs_dict(w):
    return {"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRVAL1": float(w.wcs.crval[0]), "CRVAL2": float(w.wcs.crval[1]),
            "CRPIX1": float(w.wcs.crpix[0]), "CRPIX2": float(w.wcs.crpix[1]),
            "CDELT1": float(w.wcs.cdelt[0]), "CDELT2": float(w.wcs.cdelt[1])}


def png_xy(x, y, height):
    """FITS-orientation pixel (0-based, y up) -> PNG pixel (0-based, y down)."""
    return [float(x), float(height - 1 - y)]


def limits(images, signed=False):
    """One stretch for all frames so a blink shows real change, not per-frame rescaling."""
    d = np.concatenate([i[np.isfinite(i)] for i in images]) if images else np.array([])
    if d.size == 0:
        return 0.0, 1.0
    if signed:
        return 0.0, float(np.percentile(np.abs(d), 99)) or 1.0
    lo, hi = np.percentile(d, [1, 99.5])
    return float(lo), float(hi)


def save_png(img, path, lo, hi, signed=False):
    if signed:  # diff: zero is mid-grey
        v = 0.5 + 0.5 * np.arcsinh(img / (0.1 * hi)) / np.arcsinh(10)
    else:
        v = np.arcsinh(10 * np.clip((img - lo) / (hi - lo + 1e-12), 0, 1)) / np.arcsinh(10)
    v = np.nan_to_num(np.clip(v, 0, 1), nan=0.0)
    Image.fromarray((np.flipud(v) * 255).astype(np.uint8)).save(path)


def save_flags(mask, path):
    """RGBA overlay: red where at least one quality flag is set in any exposure of the frame.
    Meaning of individual bits is not decoded here (the FITS header does not document them)."""
    rgba = np.zeros(mask.shape + (4,), np.uint8)
    rgba[mask] = (255, 0, 0, 160)
    Image.fromarray(np.flipud(rgba)).save(path)


def save_gif(paths, out, ms=600):
    ims = [Image.open(p).convert("P") for p in paths]
    ims[0].save(out, format="GIF", save_all=True, append_images=ims[1:], duration=ms, loop=0)
