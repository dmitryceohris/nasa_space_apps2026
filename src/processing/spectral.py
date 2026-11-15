"""Pairwise bandpass comparison: central wavelength and width cannot certify a match."""
import numpy as np

MIN_RESPONSE_OVERLAP = .98  # Project tolerance, not an instrument specification.


def compare_bandpasses(a,b):
    for p in (a,b):
        if any(p.get(k) is None or not np.isfinite(p[k]) or p[k] <= 0 for k in ("wave_um","bw_um")):
            return {"status":"provisional","reason":"missing_spectral_calibration","overlap":None}
    low = max(a["wave_um"]-a["bw_um"]/2,b["wave_um"]-b["bw_um"]/2)
    high = min(a["wave_um"]+a["bw_um"]/2,b["wave_um"]+b["bw_um"]/2)
    interval_overlap = max(0,high-low)/max(a["bw_um"],b["bw_um"])
    if interval_overlap < .9:
        return {"status":"rejected","reason":"wavelength_mismatch","overlap":interval_overlap}
    if not a.get("response") or not b.get("response"):
        return {"status":"provisional","reason":"spectral_response_unavailable","overlap":None,
                "interval_overlap":interval_overlap}
    responses=[]
    for p in (a,b):
        w=np.asarray(p["response"]["wave_um"],float)
        v=np.asarray(p["response"]["throughput"],float)
        if (w.ndim!=1 or v.shape!=w.shape or len(w)<3 or not np.isfinite(w).all()
                or not np.isfinite(v).all() or np.any(np.diff(w)<=0) or np.any(v<0)
                or np.trapezoid(v,w)<=0 or not p["response"].get("source")):
            return {"status":"provisional","reason":"invalid_spectral_response","overlap":None}
        responses.append((w,v/np.trapezoid(v,w)))
    grid=np.unique(np.concatenate([r[0] for r in responses]))
    values=[np.interp(grid,w,v,left=0,right=0) for w,v in responses]
    overlap=float(np.trapezoid(np.minimum(*values),grid))
    passed=overlap>=MIN_RESPONSE_OVERLAP
    return {"status":"validated" if passed else "rejected",
            "reason":None if passed else "wavelength_mismatch","overlap":overlap}


def compatible_groups(points):
    """Complete-link groups: A~B and B~C does not imply A~C."""
    groups=[]
    for p in sorted(points,key=lambda p:(p["wave_um"],p["mjd"],p["obs_id"])):
        if compare_bandpasses(p,p)["status"]!="validated":
            continue
        for group in groups:
            if all(compare_bandpasses(p,q)["status"]=="validated" for q in group):
                group.append(p)
                break
        else:
            groups.append([p])
    return groups
