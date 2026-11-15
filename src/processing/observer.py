"""Observer geometry with explicit axis-frame provenance and a bounded fallback."""
import numpy as np
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, GCRS, ICRS, get_body_barycentric
from astropy.time import Time


def observer_position(exp, *, spacecraft_frame=None):
    """Barycentric ICRS-aligned observer position in AU.

    GEOCENTER identifies an origin, not axis orientation. Never silently treat
    unknown spacecraft axes as ICRS. A caller may supply a documented frame.
    """
    t=Time(exp.mjd,format="mjd",scale="utc")
    earth=get_body_barycentric("earth",t).xyz.to_value(u.au)
    h=exp.header or {}
    keys=("X_SC","Y_SC","Z_SC")
    if not all(k in h for k in keys):
        return {"position_au":earth.tolist(),"status":"provisional","reason":"spacecraft_position_missing",
                "orbit_bound_au":None,"time_mjd_utc":exp.mjd}
    xyz=np.array([h[k] for k in keys],float)
    if not np.isfinite(xyz).all(): raise ValueError("invalid_spacecraft_position")
    frame=spacecraft_frame or h.get("XYZ_SC_FRAME")
    origin=str(h.get("XYZ_SC_SYSTEM","")).strip().upper()
    bound=float(np.linalg.norm(xyz)*u.km.to(u.au)) if origin=="GEOCENTER" else None
    if origin!="GEOCENTER" or frame not in ("GCRS","ICRS_ALIGNED"):
        return {"position_au":earth.tolist(),"status":"provisional","reason":"spacecraft_axes_unverified",
                "orbit_bound_au":bound,"time_mjd_utc":exp.mjd}
    if frame=="GCRS":
        position=GCRS(CartesianRepresentation(xyz*u.km),obstime=t).transform_to(ICRS()).cartesian.xyz.to_value(u.au)
    else:
        position=earth+xyz*u.km.to(u.au)
    return {"position_au":position.tolist(),"status":"validated","reason":None,
            "frame":frame,"orbit_bound_au":bound,"time_mjd_utc":exp.mjd}


def parallax_factors(ra_deg,dec_deg,observer_au):
    e=np.asarray(observer_au,float)
    if e.shape[-1]!=3 or not np.isfinite(e).all(): raise ValueError("invalid observer positions")
    a,d=np.radians(ra_deg),np.radians(dec_deg)
    return np.stack([e[...,0]*np.sin(a)-e[...,1]*np.cos(a),
                     e[...,0]*np.cos(a)*np.sin(d)+e[...,1]*np.sin(a)*np.sin(d)-e[...,2]*np.cos(d)],axis=-1)
