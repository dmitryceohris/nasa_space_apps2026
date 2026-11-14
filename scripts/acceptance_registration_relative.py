"""Empirical registration of real SPHEREx exposures from the local cache (no IRSA calls).

Every exposure is put on a common grid alone (its own WCS), reference stars are detected with the project
centroid, then compared (a) with the same stars in the other exposures and (b) with Gaia DR3 (VizieR).
Run from the repository root: uv run --locked python <this> <out.json>
"""
import glob
import json
import sys
import warnings
from collections import defaultdict
from pathlib import Path

sys.path[:0] = [".", "src"]
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
from astropy import units as u  # noqa: E402
from astropy.coordinates import SkyCoord  # noqa: E402
from astropy.time import Time  # noqa: E402

from Exposure import load_exposure  # noqa: E402
from processing.centroid import _detect  # noqa: E402
from processing.frames import coadd, make_grid  # noqa: E402
from processing.targets import load_targets  # noqa: E402

GROUPS = {("barnard", 3, 73): 5.0, ("epsindib", 3, 73): 5.0, ("rwaur", 2, 44): 3.5, ("rdra", 2, 44): 3.5,
          ("wise0855", 6, 44): 3.5}
T = load_targets()
files = defaultdict(list)
for p in sorted(glob.glob("cache/cut_*.fits")):
    try:
        e = load_exposure(p)
    except Exception:  # noqa: BLE001
        continue
    for (name, det, n), fov in GROUPS.items():
        if e.detector == det and max(e.shape) == n and e.contains(T[name].at(e.mjd), margin_px=3):
            files[(name, det, n)].append(e)

def gaia_query(center, radius_deg):
    """Gaia DR3 via the ESA archive (VizieR answered 404 during this check)."""
    from astroquery.gaia import Gaia
    from astropy.table import Table
    q = (f"SELECT ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE phot_g_mean_mag < 18 AND "
         f"1=CONTAINS(POINT(ra, dec), CIRCLE({center.ra.deg}, {center.dec.deg}, {radius_deg}))")
    t = Gaia.launch_job(q).get_results()
    return Table({"RA_ICRS": t["ra"], "DE_ICRS": t["dec"], "pmRA": t["pmra"], "pmDE": t["pmdec"]})


VZ = True


def robust(v):
    v = np.asarray(v, float)
    if v.size == 0:
        return None, None
    m = float(np.median(v))
    return m, (float(1.4826 * np.median(np.abs(v - m))) if v.size >= 3 else None)


out = []
for key, exps in files.items():
    name, det, n = key
    fov = GROUPS[key]
    # one exposure per obs_id
    seen, uniq = set(), []
    for e in sorted(exps, key=lambda e: e.mjd):
        if e.obs_id not in seen:
            seen.add(e.obs_id)
            uniq.append(e)
    t = T[name]
    center = t.at(np.mean([e.mjd for e in uniq]))
    wcs, shape = make_grid(center, fov)
    gaia = None
    if VZ is not None:
        try:
            gaia = gaia_query(center, fov / 60 * 0.75)
        except Exception as ex:  # noqa: BLE001
            print("gaia failed", ex, flush=True)
    rows = []
    dets = []
    for e in uniq:
        f = coadd([e], wcs, shape, t.at(e.mjd))
        if f is None:
            continue
        px, py = (float(v) for v in wcs.world_to_pixel(t.at(e.mjd)))
        objs, _, _ = _detect(f.image, 5.0, 3, 1.0, f.contam, None)
        refs = [o for o in objs if o["quality"] == "good" and o["snr"] >= 15 and (o["axis_ratio"] or 9) <= 1.4
                and np.hypot(o["xy"][0] - px, o["xy"][1] - py) > 6]
        dets.append((e, f, refs, (px, py), objs))
    if not dets:
        continue
    k0 = int(np.argmax([len(d[2]) for d in dets]))
    for e, f, refs, (px, py), objs in dets:
        dx, dy = [], []
        for o in refs:
            cand = [q for q in dets[k0][2] if np.hypot(q["xy"][0] - o["xy"][0], q["xy"][1] - o["xy"][1]) < 1.5]
            if len(cand) == 1:
                dx.append(o["xy"][0] - cand[0]["xy"][0])
                dy.append(o["xy"][1] - cand[0]["xy"][1])
        rel = dict(zip(("dx", "sx"), robust(dx))) | dict(zip(("dy", "sy"), robust(dy))) | {"n": len(dx)}
        ab = {"n": 0}
        if gaia is not None and len(gaia):
            dt = Time(e.mjd, format="mjd").jyear - 2016.0
            gde = np.asarray(gaia["DE_ICRS"], float)
            gra = np.asarray(gaia["RA_ICRS"], float) + np.nan_to_num(np.asarray(gaia["pmRA"], float)) / 3.6e6 * dt \
                / np.cos(np.radians(gde))
            gde = gde + np.nan_to_num(np.asarray(gaia["pmDE"], float)) / 3.6e6 * dt
            gx, gy = wcs.world_to_pixel(SkyCoord(gra * u.deg, gde * u.deg))
            gdx, gdy = [], []
            for o in refs:
                dd = np.hypot(gx - o["xy"][0], gy - o["xy"][1])
                srt = np.sort(dd)
                if srt[0] < 1.5 and (len(srt) < 2 or srt[1] > 3.0):
                    j = int(np.argmin(dd))
                    gdx.append(o["xy"][0] - gx[j])
                    gdy.append(o["xy"][1] - gy[j])
            ab = dict(zip(("dx", "sx"), robust(gdx))) | dict(zip(("dy", "sy"), robust(gdy))) | {"n": len(gdx)}
        # the target itself: nearest detection to its prediction (for the report only, not a pipeline result)
        near = sorted(objs, key=lambda o: np.hypot(o["xy"][0] - px, o["xy"][1] - py))
        tgt = None
        if near and np.hypot(near[0]["xy"][0] - px, near[0]["xy"][1] - py) < 3:
            o = near[0]
            tgt = {"dx": o["xy"][0] - px, "dy": o["xy"][1] - py, "quality": o["quality"], "flags": o["flags"],
                   "snr": o["snr"]}
        wl = f.wave_um
        rows.append({"obs_id": e.obs_id, "mjd": round(e.mjd, 4), "date": f.date,
                     "wave_um": None if not np.isfinite(wl) else round(wl, 4), "n_refs": len(refs),
                     "relative": rel, "gaia": ab, "target_vs_prediction": tgt})
    out.append({"target": name, "det": det, "native_size": n, "grid_fov_arcmin": fov,
                "pixel_arcsec": float(abs(wcs.wcs.cdelt[1]) * 3600), "reference_exposure": dets[k0][0].obs_id,
                "n_gaia": 0 if gaia is None else len(gaia), "exposures": rows})
    print(f"== {name} det{det} ({len(rows)} exposures, gaia {0 if gaia is None else len(gaia)})", flush=True)
    f3 = lambda v: "  -  " if v is None else f"{v:+.2f}"  # noqa: E731
    for r in rows:
        a, b, tg = r["relative"], r["gaia"], r["target_vs_prediction"]
        print(f"  {r['date']} {r['obs_id']:20s} λ={r['wave_um']} refs={r['n_refs']:2d} | rel n={a['n']:2d} "
              f"{f3(a.get('dx'))},{f3(a.get('dy'))} σ {f3(a.get('sx'))},{f3(a.get('sy'))} | gaia n={b['n']:2d} "
              f"{f3(b.get('dx'))},{f3(b.get('dy'))} σ {f3(b.get('sx'))},{f3(b.get('sy'))} | target "
              + ("-" if tg is None else f"{tg['dx']:+.2f},{tg['dy']:+.2f} {tg['quality']}"), flush=True)

Path(sys.argv[1]).write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("done", flush=True)
