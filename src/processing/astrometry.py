"""Independent native detections and reference-star registration with covariance."""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist

from .centroid import detect_all
from .psf import fit_epsf
from .validation import quality


def native_detections(exp):
    """Catalogue-free candidate detection, followed by free-position native ePSF fits."""
    image = np.where(exp.bad_mask(True), np.nan, exp.background_subtracted())
    result = detect_all(image,saturated=exp.bad_mask(True))
    fits=[]
    for candidate in result:
        if candidate["quality"] == "poor":
            continue
        p=fit_epsf(exp,exp.pixel_to_sky(*candidate["xy"]),fit_position=True)
        if p is not None:
            p["detection_xy"]=candidate["xy"]
            fits.append(p)
    return fits


def match_reference_stars(detected_xy, reference_xy, max_distance=1.5):
    """Unique reciprocal matches only. Reference positions must be propagated to the exposure."""
    measured=np.asarray(detected_xy,float).reshape(-1,2)
    reference=np.asarray(reference_xy,float).reshape(-1,2)
    if not len(measured) or not len(reference): return []
    distances=cdist(measured,reference)
    pairs=[]
    for i,row in enumerate(distances):
        close=np.flatnonzero(row<max_distance)
        if len(close)==1 and np.count_nonzero(distances[:,close[0]]<max_distance)==1:
            pairs.append((i,int(close[0])))
    return pairs


def fit_registration(measured_xy, reference_xy, covariance_px2, *, references_propagated=False):
    """Robust affine registration; parameter covariance and residual floor are propagated.

    Six coefficients; require >=8 unique stars and a nondegenerate spatial layout.
    Coordinates are in the same tangent-plane pixel convention, never raw RA/Dec.
    """
    xy,ref,cov=map(lambda v:np.asarray(v,float),(measured_xy,reference_xy,covariance_px2))
    if xy.ndim!=2 or xy.shape[1:]!=(2,) or ref.shape!=xy.shape or cov.shape!=(len(xy),2,2):
        raise ValueError("registration_shape_mismatch")
    failed={"matrix":None,"covariance":None,"residual_covariance_px2":None,
            "n_references":len(xy),"validation":quality({"astrometric_registration_failed":False},["astrometric_registration_failed"])}
    if len(xy)<8 or not np.isfinite(xy).all() or not np.isfinite(ref).all() or not np.isfinite(cov).all(): return failed
    if any(np.min(np.linalg.eigvalsh(c))<=0 for c in cov): return failed
    center=xy.mean(axis=0)
    scale=max(float(np.std(xy)),1.)
    X=np.column_stack([np.ones(len(xy)),(xy-center)/scale])
    n=len(xy)
    A=np.zeros((n,2,6)); A[:,0,:3]=X; A[:,1,3:]=X
    L=np.linalg.cholesky(cov)
    Aw=np.linalg.solve(L,A); bw=np.linalg.solve(L,ref[...,None])[...,0]
    keep=np.ones(n,bool)
    for _ in range(4):
        a,b=Aw[keep].reshape(-1,6),bw[keep].reshape(-1)
        if keep.sum()<8 or np.linalg.matrix_rank(a)<6: return failed
        coef=np.linalg.lstsq(a,b,rcond=None)[0]
        residual=ref-np.einsum('nij,j->ni',A,coef)
        delta=residual-np.median(residual[keep],axis=0)
        scatter=np.maximum(1.4826*np.median(np.abs(delta[keep]),axis=0),np.sqrt(np.median(np.diagonal(cov,axis1=1,axis2=2),axis=0)))
        revised=np.sum((delta/scatter)**2,axis=1)<=16
        if np.array_equal(revised,keep): break
        keep=revised
    if keep.sum()<8: return failed
    a,b=Aw[keep].reshape(-1,6),bw[keep].reshape(-1)
    if np.linalg.matrix_rank(a)<6: return failed
    coef=np.linalg.lstsq(a,b,rcond=None)[0]
    residual=ref-np.einsum('nij,j->ni',A,coef)
    scatter=np.cov(residual[keep].T,ddof=3)
    pcov=np.linalg.inv(a.T@a)
    # Estimate a residual floor separately; do not pretend more stars eliminate systematics.
    v=quality({"astrometric_registration_failed":True,"reference_motion_unverified":references_propagated},
              ["astrometric_registration_failed"])
    return {"matrix":coef.reshape(2,3).tolist(),"covariance":pcov.tolist(),
            "center":center.tolist(),"scale":scale,"n_references":int(keep.sum()),
            "used_indices":np.flatnonzero(keep).tolist(),"residuals_px":residual.tolist(),
            "residual_covariance_px2":scatter.tolist(),"validation":v}


def apply_registration(registration, xy, covariance_px2):
    if registration["matrix"] is None: raise ValueError("astrometric_registration_failed")
    mat=np.asarray(registration["matrix"])
    z=(np.asarray(xy)-registration["center"])/registration["scale"]
    X=np.r_[1,z]
    B=np.zeros((2,6)); B[0,:3]=X; B[1,3:]=X
    J=mat[:,1:]/registration["scale"]
    cov=J@np.asarray(covariance_px2)@J.T+B@np.asarray(registration["covariance"])@B.T
    cov+=np.asarray(registration["residual_covariance_px2"])
    return (mat@X).tolist(),cov.tolist()


def register_catalogue(exp, detections, reference_sky, reference_covariance_px2=None):
    """Register against an explicit per-epoch catalogue, excluding the target upstream.

    reference_sky is already propagated (including annual parallax when material).
    Missing catalogue covariance cannot produce validated registration.
    """
    clean=[p for p in detections if p.get("ok") and p.get("covariance_xy_px2") is not None]
    measured=[p["xy"] for p in clean]
    expected=[list(exp.sky_to_pixel(c)) for c in reference_sky]
    matches=match_reference_stars(measured,expected)
    cov=[]
    for i,j in matches:
        c=np.asarray(clean[i]["covariance_xy_px2"])
        if reference_covariance_px2 is not None: c=c+reference_covariance_px2[j]
        cov.append(c)
    if not matches:
        return fit_registration(np.empty((0,2)),np.empty((0,2)),np.empty((0,2,2)))
    return fit_registration([measured[i] for i,j in matches],[expected[j] for i,j in matches],cov,
                            references_propagated=reference_covariance_px2 is not None)
