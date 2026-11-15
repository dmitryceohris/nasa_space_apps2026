import sys, glob, json, warnings
sys.path[:0] = [".", "src"]; warnings.filterwarnings("ignore")
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.gaia import Gaia
from Exposure import load_exposure
from processing.centroid import _detect
from processing.frames import coadd, make_grid
from processing.targets import load_targets
T = load_targets()
GROUPS = (("epsindib", 3, 73, 5.0), ("barnard", 3, 73, 5.0), ("barnard", 3, 44, 3.5), ("rwaur", 2, 44, 3.5),
          ("rdra", 2, 44, 3.5), ("wise0855", 6, 44, 3.5))
allexp = [load_exposure(p) for p in sorted(glob.glob("cache/cut_*.fits"))]
out = []
for name, det, n, fov in GROUPS:
    seen, exps = set(), []
    for e in sorted(allexp, key=lambda e: e.mjd):
        if e.detector == det and max(e.shape) == n and e.contains(T[name].at(e.mjd), margin_px=3) and e.obs_id not in seen:
            seen.add(e.obs_id); exps.append(e)
    t = T[name]; center = t.at(np.mean([e.mjd for e in exps]))
    wcs, shape = make_grid(center, fov)
    g = Gaia.launch_job(f"SELECT ra, dec, pmra, pmdec FROM gaiadr3.gaia_source WHERE phot_g_mean_mag < 18 AND 1=CONTAINS(POINT(ra, dec), CIRCLE({center.ra.deg}, {center.dec.deg}, {fov/60*0.75}))").get_results()
    per, pooled = [], []
    for e in exps:
        f = coadd([e], wcs, shape, t.at(e.mjd))
        if f is None: continue
        px, py = (float(v) for v in wcs.world_to_pixel(t.at(e.mjd)))
        objs, _, _ = _detect(f.image, 5.0, 3, 1.0, f.contam, None)
        refs = [o for o in objs if o["snr"] >= 15 and not {"touches_edge", "touches_nodata"} & set(o["flags"])
                and np.hypot(o["xy"][0]-px, o["xy"][1]-py) > 6]
        dt = Time(e.mjd, format="mjd").jyear - 2016.0
        gde = np.asarray(g["dec"], float); gra = np.asarray(g["ra"], float) + np.nan_to_num(np.asarray(g["pmra"], float))/3.6e6*dt/np.cos(np.radians(gde))
        gde = gde + np.nan_to_num(np.asarray(g["pmdec"], float))/3.6e6*dt
        gx, gy = wcs.world_to_pixel(SkyCoord(gra*u.deg, gde*u.deg))
        dx, dy = [], []
        for o in refs:
            dd = np.hypot(gx - o["xy"][0], gy - o["xy"][1]); j = np.argsort(dd)
            if dd[j[0]] < 1.5 and (len(j) < 2 or dd[j[1]] > 1.5):
                dx.append(o["xy"][0]-gx[j[0]]); dy.append(o["xy"][1]-gy[j[0]])
        pooled += list(zip(dx, dy))
        if len(dx) >= 3:
            per.append({"obs_id": e.obs_id, "date": f.date, "n": len(dx), "dx": float(np.median(dx)), "dy": float(np.median(dy))})
    pooled = np.array(pooled)
    rec = {"target": name, "det": det, "native": n, "fov": fov, "n_exposures": len(exps), "n_gaia": len(g),
           "n_matches": len(pooled), "per_exposure": per}
    if len(pooled):
        rec |= {"pooled_dx_median": float(np.median(pooled[:,0])), "pooled_dy_median": float(np.median(pooled[:,1])),
                "pooled_scatter_x": float(1.4826*np.median(np.abs(pooled[:,0]-np.median(pooled[:,0])))),
                "pooled_scatter_y": float(1.4826*np.median(np.abs(pooled[:,1]-np.median(pooled[:,1]))))}
    if per:
        off = np.array([np.hypot(p["dx"], p["dy"]) for p in per])
        rec |= {"exposure_offset_median": float(np.median(off)), "exposure_offset_p90": float(np.percentile(off, 90)),
                "exposure_offset_max": float(off.max())}
    out.append(rec)
    print(f"{name:9s} det{det} n{n}: exp {len(exps):2d}, gaia {len(g):3d}, matches {len(pooled):4d}" +
          (f", pooled median ({rec['pooled_dx_median']:+.3f},{rec['pooled_dy_median']:+.3f}) px, star scatter ({rec['pooled_scatter_x']:.3f},{rec['pooled_scatter_y']:.3f}) px" if len(pooled) else "") +
          (f"; {len(per)} exp with >=3: |offset| median {rec['exposure_offset_median']:.3f} p90 {rec['exposure_offset_p90']:.3f} max {rec['exposure_offset_max']:.3f} px" if per else "; no exposure with >=3 matches"), flush=True)
json.dump(out, open(sys.argv[1], "w", encoding="utf-8"), ensure_ascii=False, indent=1)
