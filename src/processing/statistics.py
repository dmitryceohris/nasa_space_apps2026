"""Covariance-aware constant flux test; correlations are explicit inputs."""
import numpy as np
from scipy.stats import chi2


def constant_flux(points, *, covariance_ujy2=None, calibration_correlation=None):
    from .validation import accepted
    from .spectral import compare_bandpasses
    if len(points)<2 or not all(accepted(p) for p in points):
        return {"status":"provisional","reason":"insufficient_validated_photometry"}
    identities=[(p['obs_id'],p.get('det')) for p in points]
    if len(set(identities))!=len(identities):
        return {"status":"rejected","reason":"duplicate_physical_observation"}
    if any(compare_bandpasses(a,b)["status"]!="validated" for i,a in enumerate(points) for b in points[i+1:]):
        return {"status":"rejected","reason":"wavelength_mismatch"}
    f=np.array([p["flux_ujy"] for p in points],float)
    n=len(f)
    if covariance_ujy2 is None:
        if calibration_correlation is None:
            return {"status":"provisional","reason":"calibration_correlation_unknown"}
        if not 0<=calibration_correlation<=1: raise ValueError("calibration correlation must be in [0,1]")
        stat=np.array([p["stat_err_ujy"] for p in points])
        sys=np.array([p["systematic_err_ujy"] for p in points])
        C=np.diag(stat**2+(1-calibration_correlation)*sys**2)+calibration_correlation*np.outer(sys,sys)
    else: C=np.asarray(covariance_ujy2,float)
    if C.shape!=(n,n) or not np.isfinite(C).all() or not np.allclose(C,C.T): raise ValueError("invalid flux covariance")
    np.linalg.cholesky(C)
    one=np.ones(n); invone=np.linalg.solve(C,one)
    mean=float(invone@f/(one@invone))
    residual=f-mean
    chisq=float(residual@np.linalg.solve(C,residual))
    return {"status":"validated","weighted_mean_ujy":mean,"mean_err_ujy":float(np.sqrt(1/(one@invone))),
            "chi2":chisq,"dof":n-1,"p_value":float(chi2.sf(chisq,n-1)),"covariance_ujy2":C.tolist()}
