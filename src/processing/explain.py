"""What happens inside a rectangle of the field: facts for the "Explain this region" card (docs/EXPLAIN_TEXTS.md).

Input: the per-epoch grids saved next to an /api/view result (api/jobs._save_arrays). Output: plain facts
(verdict, numbers, caveat conditions); the wording lives in api/explain_text.py. Same input -> same output.

Verdicts, checked in this order (the first that fits wins):
  not_enough_data  fewer than 2 epochs with data in the area (one card for the whole area)
  cant_tell        saturated core (checked first: a saturated star is real, its bloom flags are not an artefact)
  artefact         every sign of change sits on contamination-flagged pixels
  cant_tell        two sources too close / source touching the edge of the data
  moving           a +/- residual pair in >= 2 epochs and a significant linear shift of the centre
  seen_once        a residual in exactly one epoch, nothing at that place in the others
  brightness       a significant flux change at a fixed place
  nothing          no change above the noise (one card for the whole area, only when no other card exists)
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from .centroid import detect
from .frames import FLAG_CODE_NAMES, difference_details

MIN_COVER = 0.5  # an epoch has data in the area when at least this fraction of its pixels is valid
FULL_COVER = 0.9  # below this an epoch counts as partially covered (caveat)
MATCH_PX = 3.0  # detections closer than this in different epochs belong to the same source
PAIR_PX = 4.0  # + and - residuals this close to a source form its pair
AMBIGUOUS_PX = 3.0  # two sources closer than this: we cannot tell which one changes
AP_PX = 2.0  # radius of the disc used for flags and fluxes
REGISTRATION_PX = 0.4  # added to every centroid error: largest WCS offset found against Gaia DR3
MOTION_SIGMA = 3.0  # a shift is real when it exceeds this many errors
FLUX_SIGMA = 5.0  # a brightness change is real when it exceeds this many errors
ARTEFACT_FRACTION = 0.5  # share of contamination-flagged pixels under a residual that makes it an artefact
PREDICTION_NEAR_PX = 5.0  # a source this close to the predicted path is compared with the prediction
PREDICTION_AGREE_PX = 1.0  # same gate as /api/view (api/jobs.ASSOC_PX)
FAST_ARCSEC_YR = 0.5
SMALL_SHIFT_PX = 2.0
PARALLAX_SHARE = 0.2

SATURATION = ("OVERFLOW", "NONLINEAR")
SAT_BITS = sum(1 << FLAG_CODE_NAMES.index(n) for n in SATURATION)
CONTAM_BITS = sum(1 << k for k, n in enumerate(FLAG_CODE_NAMES) if n not in SATURATION)
DIRECTIONS = ("north", "north-east", "east", "south-east", "south", "south-west", "west", "north-west")


def _disc(shape, xy, r=AP_PX):
    yy, xx = np.indices(shape)
    return np.hypot(xx - xy[0], yy - xy[1]) <= r


def _near(objs, xy, r):
    return [o for o in objs if np.hypot(o["xy"][0] - xy[0], o["xy"][1] - xy[1]) <= r]


def _nearest(objs, xy, r):
    near = sorted(_near(objs, xy, r), key=lambda o: np.hypot(o["xy"][0] - xy[0], o["xy"][1] - xy[1]))
    return near[0] if near else None


def _flag_names(codes):
    names = []
    for k, n in enumerate(FLAG_CODE_NAMES):
        if n not in SATURATION and np.any(codes & (1 << k)):
            names.append(n)
    return names


def _direction(east, north):
    angle = np.degrees(np.arctan2(east, north)) % 360  # 0 = north, 90 = east
    return DIRECTIONS[int(((angle + 22.5) % 360) // 45)]


def _fit(t, v, s):
    """Weighted straight line v = a + b t. -> (a, b, sigma_b)."""
    w = 1 / s ** 2
    tm = np.sum(w * t) / np.sum(w)
    sxx = np.sum(w * (t - tm) ** 2)
    b = np.sum(w * (t - tm) * v) / sxx
    a = np.sum(w * v) / np.sum(w) - b * tm
    return a, b, 1 / np.sqrt(sxx)


class _Region:
    def __init__(self, arrays, box, pixel_scale_arcsec, target, parallax_mas):
        self.img, self.var = arrays["image"].astype(float), arrays["variance"].astype(float)
        self.codes = arrays["codes"].astype(np.int64)
        self.mjd = np.asarray(arrays["mjd"], float)
        self.wave_min, self.wave_max = np.asarray(arrays["wave_min"], float), np.asarray(arrays["wave_max"], float)
        self.channel = float(arrays["channel_width"])
        self.pred = np.asarray(arrays["predicted"], float)
        self.scale, self.target, self.parallax = pixel_scale_arcsec, target, parallax_mas
        n, h, w = self.img.shape
        x0, y0, x1, y1 = box
        self.x0, self.x1 = max(0, int(np.floor(x0))), min(w - 1, int(np.ceil(x1)))
        self.y0, self.y1 = max(0, int(np.floor(y0))), min(h - 1, int(np.ceil(y1)))
        self.shape = (h, w)

    def inside(self, xy):
        return self.x0 - 0.5 <= xy[0] <= self.x1 + 0.5 and self.y0 - 0.5 <= xy[1] <= self.y1 + 0.5

    def cut(self, a):
        return a[..., self.y0:self.y1 + 1, self.x0:self.x1 + 1]


def analyse(arrays: dict, box, pixel_scale_arcsec: float, target: str | None = None,
            parallax_mas: float | None = None) -> dict:
    """box: x0, y0, x1, y1 in grid pixels (0-based, y up, inclusive). -> {n_epochs, n_epochs_with_data,
    n_partial, items: [{verdict, xy, ...facts}]}; see the module docstring."""
    R = _Region(arrays, box, pixel_scale_arcsec, target, parallax_mas)
    n = R.img.shape[0]
    # saturated pixels hold no value but are not missing data: they count as covered
    has = np.isfinite(R.img) | ((R.codes & SAT_BITS) > 0)
    cover = R.cut(has).reshape(n, -1).mean(axis=1) if R.x1 >= R.x0 and R.y1 >= R.y0 else np.zeros(n)
    use = np.flatnonzero(cover >= MIN_COVER)
    base = {"n_epochs": int(n), "n_epochs_with_data": int(use.size), "n_partial": int(np.sum(cover < FULL_COVER))}
    if use.size < 2:
        return {**base, "items": [{"verdict": "not_enough_data", "n_epochs_with_data": int(use.size)}]}

    img, var, codes, mjd = R.img[use], R.var[use], R.codes[use], R.mjd[use]
    dd = difference_details([SimpleNamespace(image=img[j], variance=var[j]) for j in range(use.size)])
    diffs, template = dd["diffs"], dd["template"]
    sat = (codes & SAT_BITS) > 0
    epoch_objs = [detect(img[j], saturated=sat[j]) for j in range(use.size)]
    diff_objs = [detect(diffs[j], signed=True, saturated=sat[j]) for j in range(use.size)]
    template_objs = [o for o in detect(template, saturated=sat.any(axis=0)) if R.inside(o["xy"])]

    def contam_share(j, xy):
        d = _disc(R.shape, xy)
        return float(np.mean((codes[j][d] & CONTAM_BITS) > 0)), codes[j][d]

    def flagged(j, xy):
        return bool(np.any(codes[j][_disc(R.shape, xy)]))

    def saturated(j, xy):
        return bool(np.any(codes[j][_disc(R.shape, xy)] & SAT_BITS))

    def evidence(xy):
        return [(j, o) for j in range(use.size) for o in _near(diff_objs[j], xy, PAIR_PX)]

    def artefact(ev):
        shares = [contam_share(j, o["xy"]) for j, o in ev]
        if ev and all(s >= ARTEFACT_FRACTION for s, _ in shares):
            return _flag_names(np.bitwise_or.reduce(np.concatenate([c for _, c in shares])))
        return None

    items, taken = [], set()
    # residuals that belong to no persistent source: transients (one epoch) or sources missing in the template
    groups = []
    for j in range(use.size):
        for o in diff_objs[j]:
            if o["sign"] != 1 or not R.inside(o["xy"]) or _near(template_objs, o["xy"], MATCH_PX):
                continue
            g = next((g for g in groups if np.hypot(*(np.subtract(g["xy"], o["xy"]))) <= MATCH_PX), None)
            if g is None:
                groups.append({"xy": list(o["xy"]), "hits": [(j, o)]})
            elif j not in {h[0] for h in g["hits"]}:
                g["hits"].append((j, o))
    sources = [{"xy": o["xy"], "obj": o} for o in template_objs]
    for g in groups:
        if len(g["hits"]) == 1:
            j, o = g["hits"][0]
            names = artefact([(j, o)])
            if names is not None:
                items.append({"verdict": "artefact", "xy": o["xy"], "flag_names": names})
            else:
                items.append({"verdict": "seen_once", "xy": o["xy"], "mjd_seen": float(mjd[j]),
                              "n_other": int(use.size - 1)})
        else:
            sources.append({"xy": list(np.mean([o["xy"] for _, o in g["hits"]], axis=0)), "obj": None})

    for si, s in enumerate(sources):
        if si in taken:
            continue
        xy, obj = s["xy"], s["obj"]
        ev = evidence(xy)
        # A saturated source is real: its bloom/halo flags must not turn it into an "artefact".
        n_sat = sum(saturated(j, xy) for j in range(use.size))
        if (obj is not None and "saturated_core" in obj["flags"]) or n_sat > use.size / 2:
            items.append({"verdict": "cant_tell", "reason": "saturated", "xy": xy})
            continue
        names = artefact(ev)
        if names is not None:
            items.append({"verdict": "artefact", "xy": xy, "flag_names": names})
            continue
        if not ev:
            continue
        close = [k for k, t in enumerate(sources) if k != si and k not in taken
                 and np.hypot(t["xy"][0] - xy[0], t["xy"][1] - xy[1]) < AMBIGUOUS_PX]
        if close:
            taken.update(close)
            items.append({"verdict": "cant_tell", "reason": "ambiguous", "xy": xy})
            continue
        if obj is not None and {"touches_edge", "touches_nodata"} & set(obj["flags"]):
            items.append({"verdict": "cant_tell", "reason": "edge", "xy": xy})
            continue
        item = _moving(R, use, img, mjd, epoch_objs, diff_objs, xy, flagged) or \
            _brightness(R, use, img, var, mjd, xy, flagged)
        if item is not None:
            items.append(item)

    if not items:
        items = [{"verdict": "nothing"}]
    return {**base, "items": items}


def _wavelengths(R, use, epochs):
    lo, hi = R.wave_min[use][epochs], R.wave_max[use][epochs]
    if not np.isfinite(lo).any():
        return {}
    w0, w1 = float(np.nanmin(lo)), float(np.nanmax(hi))
    width = R.channel if np.isfinite(R.channel) and R.channel > 0 else 0.05 * (w0 + w1) / 2
    return {"wave_min": w0, "wave_max": w1, "different_wavelengths": bool(w1 - w0 > width)}


def _moving(R, use, img, mjd, epoch_objs, diff_objs, xy, flagged):
    rows, n_flagged = [], 0
    for j in range(use.size):
        o = _nearest(epoch_objs[j], xy, MATCH_PX)
        if o is None:
            continue
        if flagged(j, o["xy"]):
            n_flagged += 1
            continue
        if o["quality"] == "poor":
            continue
        rows.append((j, o))
    n_pair = sum(1 for j in range(use.size)
                 if any(o["sign"] == 1 for o in _near(diff_objs[j], xy, PAIR_PX))
                 and any(o["sign"] == -1 for o in _near(diff_objs[j], xy, PAIR_PX)))
    if len(rows) < 2 or n_pair < 2:
        return None
    t = np.array([(mjd[j] - mjd[rows[0][0]]) / 365.25 for j, _ in rows])
    if np.ptp(t) <= 0:
        return None
    x = np.array([o["xy"][0] for _, o in rows])
    y = np.array([o["xy"][1] for _, o in rows])
    sx = np.array([np.hypot(o["err_px"][0], REGISTRATION_PX) for _, o in rows])
    sy = np.array([np.hypot(o["err_px"][1], REGISTRATION_PX) for _, o in rows])
    ax, vx, svx = _fit(t, x, sx)
    ay, vy, svy = _fit(t, y, sy)
    pm = float(np.hypot(vx, vy))
    pm_err = float(np.sqrt((vx * svx) ** 2 + (vy * svy) ** 2) / pm) if pm > 0 else float("inf")
    if not pm > MOTION_SIGMA * pm_err:
        return None
    span = float(np.ptp(t))
    shift_px = pm * span
    item = {"verdict": "moving", "xy": xy, "n_pair": int(n_pair), "n_epochs": int(use.size),
            "mjd_first": float(mjd[rows[0][0]]), "mjd_last": float(mjd[rows[-1][0]]),
            "shift_px": shift_px, "shift_arcsec": shift_px * R.scale,
            "pm_arcsec_yr": pm * R.scale, "pm_err_arcsec_yr": pm_err * R.scale,
            "direction": _direction(-vx, vy),  # grid: x grows to the west, y to the north
            "n_flagged": int(n_flagged), "small_shift": bool(shift_px < SMALL_SHIFT_PX)}
    item.update(_wavelengths(R, use, [j for j, _ in rows]))

    pred = R.pred[use]
    have = np.isfinite(pred).all(axis=1)
    item["prediction"] = "none"
    if R.target and have.any():
        tt = (mjd[have] - mjd[rows[0][0]]) / 365.25
        fit = np.column_stack([ax + vx * tt, ay + vy * tt])
        if np.min(np.hypot(*(pred[have] - np.asarray(xy)).T)) <= PREDICTION_NEAR_PX:
            offset = float(np.max(np.hypot(*(fit - pred[have]).T)))
            item.update(prediction="agrees" if offset <= PREDICTION_AGREE_PX else "disagrees",
                        offset_px=offset, target=R.target)
            if R.parallax and R.parallax > 0:
                item["parallax_noticeable"] = bool(R.parallax / 1000 >= PARALLAX_SHARE * item["shift_arcsec"])
    if item["prediction"] == "none":
        item["fast"] = bool(item["pm_arcsec_yr"] >= FAST_ARCSEC_YR)
    return item


def _brightness(R, use, img, var, mjd, xy, flagged):
    disc = _disc(R.shape, xy)
    rows, n_flagged = [], 0
    for j in range(use.size):
        if not np.isfinite(img[j][disc]).all():
            continue
        if flagged(j, xy):
            n_flagged += 1
            continue
        rows.append((j, float(np.sum(img[j][disc])), float(np.sum(var[j][disc]))))
    if len(rows) < 2:
        return None
    jb, fb, vb = max(rows, key=lambda r: r[1])
    jf, ff, vf = min(rows, key=lambda r: r[1])
    if fb <= 0 or not (fb - ff) > FLUX_SIGMA * np.sqrt(vb + vf):
        return None
    item = {"verdict": "brightness", "xy": xy, "change_pct": 100 * (fb - ff) / fb,
            "mjd_bright": float(mjd[jb]), "mjd_faint": float(mjd[jf]), "n_flagged": int(n_flagged)}
    item.update(_wavelengths(R, use, [j for j, _, _ in rows]))
    return item
