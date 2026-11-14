"""Object centre in a rendered frame, with no labels and no assumption about size.

Steps: robust background and noise (median / MAD) -> light smoothing -> threshold -> connected components ->
for each component a flux-weighted centroid. Saturated or masked pixels (NaN) inside an object are treated as part
of it: they count towards its area and get its brightest weight, so a saturated star is not biased by its hole.

Two candidates are reported when they differ: the LARGEST object and the CENTRAL one (the one under the image
centre, else the nearest). If they are the same object there is one candidate. Nothing found -> empty list.
Limits: the background estimate fails if one object fills most of the frame; touching objects merge.
"""
import numpy as np
from scipy import ndimage as ndi

K_SIGMA = 3.0  # detection threshold in noise sigmas; raise it for noisy frames, lower it for faint objects
MIN_AREA = 3  # smallest accepted object, pixels above threshold
SMOOTH = 1.0  # gaussian sigma in pixels applied before thresholding


def find(img, signed=False, k=K_SIGMA, min_area=MIN_AREA, smooth=SMOOTH):
    """-> list of {kind: 'largest'|'central'|'largest+central', xy: [x, y] (0-based, y up like numpy), area, flux}.
    signed=True is for diff frames: both bright and dark residuals count (|value|)."""
    a = np.abs(img) if signed else np.asarray(img, float)
    fin = np.isfinite(a)
    if fin.sum() < 16:
        return []
    f = np.where(fin, a, np.median(a[fin]))
    s = ndi.gaussian_filter(f, smooth) if smooth > 0 else f
    bg = np.median(s[fin])
    sig = 1.4826 * np.median(np.abs(s[fin] - bg))
    if not sig > 0:
        return []
    above = (s - bg) > k * sig
    lab, n = ndi.label(ndi.binary_fill_holes(above | ~fin))
    w = np.where(fin, np.clip(f - bg, 0, None), 0.0)
    yy, xx = np.indices(a.shape)
    objs = []
    for i in range(1, n + 1):
        m = lab == i
        if (m & above).sum() < min_area:  # drops lone bad pixels
            continue
        wi = np.where(m, w, 0.0)
        wi[m & ~fin] = wi.max()  # saturated / masked pixels inside the object: brightest weight
        t = wi.sum()
        if t <= 0:
            continue
        objs.append(dict(id=i, area=int(m.sum()), flux=float(t),
                         xy=[float((wi * xx).sum() / t), float((wi * yy).sum() / t)]))
    if not objs:
        return []
    ny, nx = a.shape
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    largest = max(objs, key=lambda o: (o["area"], o["flux"]))
    under = lab[int(round(cy)), int(round(cx))]
    central = next((o for o in objs if o["id"] == under), None) or min(
        objs, key=lambda o: (o["xy"][0] - cx) ** 2 + (o["xy"][1] - cy) ** 2)

    def out(o, kind):
        return {"kind": kind, "xy": o["xy"], "area": o["area"], "flux": o["flux"]}

    if largest["id"] == central["id"]:
        return [out(largest, "largest+central")]
    return [out(largest, "largest"), out(central, "central")]
