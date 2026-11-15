"""Object centre in a rendered frame, with no labels and no assumption about size.

Method (``METHOD`` below, also returned by the API so the user sees how the number was obtained):
  1. robust background and noise: median / MAD with 3-sigma clipping;
  2. masked pixels (NaN, inf) are filled with the value of the NEAREST valid pixel, only to smooth and to weigh;
  3. light gaussian smoothing -> threshold k*sigma on VALID pixels only -> holes filled -> connected components;
  4. per component a flux-weighted first moment (centroid) and its noise error, plus quality flags.

NaN handling: a masked region becomes part of an object only when the object surrounds it (a saturated core is a
hole inside the bright ring). Masked areas at the frame edge or next to an object are never objects and get no
weight; an object that touches them is flagged ``touches_nodata`` because its centre is then biased.
Masked pixels inside an object take the value of their nearest valid neighbour, so a hole does not pull the
centre (a symmetric saturated core stays symmetric, a bad pixel on a wing is not given the peak value).

Difference frames (``signed=True``): positive and negative residuals are detected SEPARATELY. A moving object
leaves a "+" where it is now and a "-" where it was; measured together with |value| they merge into one blob
whose centre lies between the two positions and means nothing. The measured position is taken from "+" only.

Choice between objects: the LARGEST and the CENTRAL one (under the image centre, else the nearest) are both
reported. When they differ, or two objects are equally close to the centre, the result is "ambiguous" and no
single position is given: the server does not guess.

Not done here (later stage): PSF fitting, deblending of touching objects, neural-network detection.
"""
import numpy as np
from scipy import ndimage as ndi

K_SIGMA = 3.0  # detection threshold in noise sigmas; raise it for noisy frames, lower it for faint objects
MIN_AREA = 3  # smallest accepted object, valid pixels above threshold
SMOOTH = 1.0  # gaussian sigma in pixels applied before thresholding
HOT_PIXEL_RATIO = 0.1  # second-brightest / brightest pixel below this (with a bright peak) -> lone hot pixel
HOT_PIXEL_SNR = 10.0  # ... only when the brightest pixel is at least this many sigma
ELONGATED_AXIS_RATIO = 1.5  # second-moment axis ratio above which an object is flagged (blend, trail, artefact)
ELONGATED_MIN_SNR = 10.0  # the shape of fainter objects is noise, they are not tested
LOW_SNR = 5.0  # objects below this S/N are listed nowhere and never chosen: noise gives such blobs at 3 sigma
TIE_PX = 0.5  # two objects whose distances to the centre differ by less than this are equally "central"
DIPOLE_OVERLAP_PX = 5.0  # in a difference frame a "-" part closer than this eats the facing wing of the "+":
# a true 1.3 px shift gave "+" and "-" centres 3.0 px apart, each ~0.9 px off; at 5.9 px the bias is < 0.4 px

METHOD = {
    "name": "threshold_moment_centroid",
    "version": 2,
    "description": "Порог k·σ над робастным фоном (медиана/MAD) после сглаживания, связные области, "
                   "центроид как взвешенный по потоку первый момент. В разностных кадрах «+» и «−» "
                   "измеряются отдельно, положение берётся по «+».",
    "params": {"k_sigma": K_SIGMA, "min_area_px": MIN_AREA, "smooth_sigma_px": SMOOTH},
    "error_model": "err_px: шум центроида из шума фона (σ·√Σ(x−x̄)²/F); систематику (PSF, блендинг, "
                   "края маски) не включает",
    "quality": {"good": "без флагов", "fair": "есть маскированные/насыщенные пиксели или вытянутая форма",
                "poor": "касается края кадра или маски"},
    "min_snr": LOW_SNR,
    "limitations": [
        "фон оценивается по всему кадру: если объект занимает большую часть поля, порог завышен",
        "соприкасающиеся объекты сливаются в один (флаг elongated), разделения нет",
        "у объекта, касающегося маски или края кадра, центр смещён внутрь видимой части",
        "одиночный горячий пиксель отбрасывается только если он ярче ~10σ и соседи почти пустые",
        "в разностном кадре «+» ближе 5 пикс. к «−» смещён от него (при сдвиге 1.3 пикс. на 0.9 пикс.): "
        "такой центр помечается overlaps_negative_residual и не считается положением",
        "объекты с S/N < 5 не выбираются (на чистом шуме порог 3σ даёт такие пятна), только считаются",
        "погрешность учитывает только шум, а не форму PSF и пиксельную дискретизацию",
    ],
    "not_implemented": ["PSF-фиттинг", "деблендинг", "нейросетевое обнаружение и центроиды"],
}


def _robust(v):
    """-> (median, sigma) of finite values with two passes of 3-sigma clipping."""
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan
    med = np.median(v)
    sig = 1.4826 * np.median(np.abs(v - med))
    for _ in range(2):
        if not sig > 0:
            break
        keep = np.abs(v - med) < 3 * sig
        if keep.sum() < 8:
            break
        med = np.median(v[keep])
        sig = 1.4826 * np.median(np.abs(v[keep] - med))
    return float(med), float(sig)


def _nearest_fill(a, fin):
    if fin.all():
        return a
    idx = ndi.distance_transform_edt(~fin, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def _detect(a, k, min_area, smooth, saturated, sign):
    """Objects of one sign in `a` (already multiplied by the sign). -> (objects, n_rejected_hot, label image)."""
    fin = np.isfinite(a)
    if fin.sum() < 16:
        return [], 0, None
    f = _nearest_fill(a, fin)
    s = ndi.gaussian_filter(f, smooth) if smooth > 0 else f
    bg, sig = _robust(np.where(fin, s, np.nan))
    bg_raw, sig_raw = _robust(np.where(fin, a, np.nan))
    if not (sig > 0 and sig_raw > 0):
        return [], 0, None
    above = fin & ((s - bg) > k * sig)
    lab, n = ndi.label(ndi.binary_fill_holes(above))  # holes = masked cores / dips enclosed by the object
    w_all = np.clip(f - bg_raw, 0, None)
    yy, xx = np.indices(a.shape)
    ny, nx = a.shape
    objs, hot = [], 0
    for i in range(1, n + 1):
        m = lab == i
        if (m & above).sum() < min_area:
            continue
        good = m & fin
        vals = np.sort(w_all[good])[::-1]
        masked = m & ~fin
        if (not masked.any() and vals.size >= 2 and vals[0] > HOT_PIXEL_SNR * sig_raw
                and vals[1] < HOT_PIXEL_RATIO * vals[0]):
            hot += 1  # one bright pixel with empty neighbours: cosmic ray / hot pixel, not an object
            continue
        w = np.where(m, w_all, 0.0)
        t = w.sum()
        if t <= 0:
            continue
        x, y = (w * xx).sum() / t, (w * yy).sum() / t
        # noise error of the first moment: every valid pixel of the object carries sigma_raw
        ex = sig_raw * np.sqrt(((xx[good] - x) ** 2).sum()) / t
        ey = sig_raw * np.sqrt(((yy[good] - y) ** 2).sum()) / t
        snr = t / (sig_raw * np.sqrt(max(good.sum(), 1)))
        # second moments -> axis ratio
        cxx = (w * (xx - x) ** 2).sum() / t
        cyy = (w * (yy - y) ** 2).sum() / t
        cxy = (w * (xx - x) * (yy - y)).sum() / t
        ev = np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]])
        axis_ratio = float(np.sqrt(ev[1] / ev[0])) if ev[0] > 1e-9 else float("inf")

        flags = []
        if masked.any():
            flags.append("masked_pixels")
            if saturated is not None and (saturated & masked).any():
                flags.append("saturated_core")
        if (ndi.binary_dilation(m) & ~fin & ~m).any():
            flags.append("touches_nodata")
        if m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any():
            flags.append("touches_edge")
        if snr >= ELONGATED_MIN_SNR and axis_ratio > ELONGATED_AXIS_RATIO:
            flags.append("elongated")
        poor = {"touches_nodata", "touches_edge"}
        quality = "poor" if poor & set(flags) else "fair" if flags else "good"
        objs.append(dict(id=i, sign=sign, area=int(m.sum()), n_masked=int(masked.sum()), flux=float(t),
                         xy=[float(x), float(y)], err_px=[float(ex), float(ey)], snr=float(snr),
                         axis_ratio=round(axis_ratio, 3) if np.isfinite(axis_ratio) else None,
                         flags=flags, quality=quality))
    return objs, hot, lab


def _split_faint(objs):
    return [o for o in objs if o["snr"] >= LOW_SNR], sum(1 for o in objs if o["snr"] < LOW_SNR)


def _choose(objs, lab, shape):
    """-> (candidates with kind, ambiguity reason or None). Same rules as before: largest and central."""
    if not objs:
        return [], None
    ny, nx = shape
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    largest = max(objs, key=lambda o: (o["area"], o["flux"]))
    under = lab[int(round(cy)), int(round(cx))]
    central = next((o for o in objs if o["id"] == under), None)
    tie = False
    if central is None:
        dist = sorted(((np.hypot(o["xy"][0] - cx, o["xy"][1] - cy), o) for o in objs), key=lambda t: t[0])
        central = dist[0][1]
        tie = len(dist) > 1 and dist[1][0] - dist[0][0] < TIE_PX
    if largest["id"] == central["id"]:
        cands = [dict(largest, kind="largest+central")]
    else:
        cands = [dict(largest, kind="largest"), dict(central, kind="central")]
    reason = None
    if tie:
        reason = "два объекта на одинаковом расстоянии от центра"
    elif len(cands) == 2:
        reason = "самый большой и центральный объекты различаются"
    return cands, reason


def _public(o):
    keys = ("kind", "sign", "xy", "err_px", "area", "n_masked", "flux", "snr", "axis_ratio", "quality", "flags")
    return {k: o[k] for k in keys}


def measure(img, signed=False, saturated=None, k=K_SIGMA, min_area=MIN_AREA, smooth=SMOOTH):
    """Full measurement of one frame.

    -> {method, status: 'ok'|'ambiguous'|'none', reason, xy, err_px, quality, flags, candidates,
        n_rejected_hot, n_low_snr}
    xy / err_px / quality are set only for status 'ok'. Coordinates are 0-based [x, y] with y up (numpy).
    signed=True: difference frame in sigma units; candidates carry sign +1 / -1, the position comes from "+".
    saturated: optional bool mask (e.g. Frame.contam) used only to label masked cores as saturated.
    """
    a = np.asarray(img, float)
    sat = None if saturated is None else np.asarray(saturated, bool)
    signs = (1, -1) if signed else (1,)
    per_sign, hot, faint, found = {}, 0, {}, {}
    for sg in signs:
        objs, h, lab = _detect(sg * a, k, min_area, smooth, sat, sg if signed else None)
        objs, faint[sg] = _split_faint(objs)
        found[sg] = objs
        hot += h
        per_sign[sg] = _choose(objs, lab, a.shape)
    pos, reason = per_sign[1]
    cands = [c for sg in signs for c in per_sign[sg][0]]
    out = {"method": METHOD["name"], "status": "none", "reason": None, "xy": None, "err_px": None,
           "quality": None, "flags": [], "candidates": [_public(c) for c in cands], "n_rejected_hot": hot,
           "n_low_snr": sum(faint.values())}
    if not pos:
        if signed and per_sign[-1][0]:
            out["reason"] = "только отрицательный остаток: объект ушёл отсюда или потускнел"
        elif faint[1]:
            out["reason"] = f"есть только слабые пятна с S/N < {LOW_SNR:g} — их не отличить от шума"
        else:
            out["reason"] = f"нет объекта выше {k:g}σ"
        return out
    if reason:
        out.update(status="ambiguous", reason=reason)
        return out
    c = pos[0]
    flags, quality = list(c["flags"]), c["quality"]
    if signed and any(np.hypot(n["xy"][0] - c["xy"][0], n["xy"][1] - c["xy"][1]) < DIPOLE_OVERLAP_PX
                      for n in found[-1]):  # every '-' object, not only the listed largest/central ones
        flags.append("overlaps_negative_residual")
        quality = "poor"
    out.update(status="ok", xy=c["xy"], err_px=c["err_px"], quality=quality, flags=flags)
    return out


def find(img, signed=False, k=K_SIGMA, min_area=MIN_AREA, smooth=SMOOTH, saturated=None):
    """-> list of candidates {kind: 'largest'|'central'|'largest+central', sign, xy, err_px, area, flux, snr,
    quality, flags, ...}. signed=True: "+" and "-" residuals are separate objects (see module docstring)."""
    return measure(img, signed, saturated, k, min_area, smooth)["candidates"]


def detect_all(img, *, saturated=None):
    """All independent positive detections, without the display's largest/central selection."""
    objects, _, _ = _detect(np.asarray(img,float),K_SIGMA,MIN_AREA,SMOOTH,saturated,None)
    objects, _ = _split_faint(objects)
    return objects
