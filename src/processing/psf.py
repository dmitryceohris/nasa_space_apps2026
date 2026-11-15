"""Native-pixel position-dependent ePSF fitting, with explicit diagnostic failures.

SPHEREx ePSFs already include the pixel response: evaluate at pixel centres, do
not convolve/integrate the pixel response a second time (IRSA ePSF tutorial §7.2).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.stats import chi2

from Exposure.calibration import parent_pixel
from .validation import photometry_quality

PSF_METHOD = {
    "name": "native_epsf_fit", "params": {"radius_px": 5, "systematic_fraction": 0.03},
    "description": "Position-dependent pixel-convolved ePSF on native pixels; variance-weighted source flux and local background plane.",
    "error_components": ["pixel_variance", "fitted_background", "flux_calibration_project_floor"],
    "limitations": ["3% calibration floor is a project assumption; inter-epoch correlations require explicit covariance",
                    "unmodelled blends and unavailable calibrations prevent validation"],
    "quality": "Validated only with calibrated wavelength, solid angle, verified time and FLAGS, S/N>=5 and passing fit diagnostics",
    "not_implemented": [],
}


def select_epsf(exp, x, y):
    if exp.epsf_data is None or exp.epsf_header is None:
        raise ValueError("missing_epsf")
    px, py = parent_pixel(exp, x, y)
    table = exp.epsf_data
    j = int(np.argmin((table["XCENTER"] - px) ** 2 + (table["YCENTER"] - py) ** 2))
    row = table[j]
    if row["NSTAR"] <= 0:
        raise ValueError("unverified_epsf_zone")
    model = np.asarray(row["EPSF"], float)
    ox, oy = int(exp.epsf_header.get("OVSMPX", 0)), int(exp.epsf_header.get("OVSMPY", 0))
    if model.ndim != 2 or not np.isfinite(model).all() or ox < 1 or oy < 1 or model.sum() <= 0:
        raise ValueError("invalid_epsf")
    return model, (ox, oy), j


def render_epsf(model, oversampling, xx, yy, x0, y0):
    """Normalize on the FULL native support before cropping or masking pixels."""
    ox, oy = oversampling
    cy, cx = (np.array(model.shape) - 1) / 2
    def evaluate(x, y):
        coords = np.array([(np.asarray(y) - y0) * oy + cy, (np.asarray(x) - x0) * ox + cx])
        return map_coordinates(model, coords, order=3, mode="constant", cval=0, prefilter=True)
    ry, rx = np.ceil(np.array(model.shape) / [oy, ox] / 2).astype(int) + 1
    fy, fx = np.mgrid[-ry:ry+1, -rx:rx+1]
    norm = evaluate(fx + round(x0), fy + round(y0)).sum()
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("invalid_epsf_normalization")
    return evaluate(xx, yy) / norm


def fit_epsf(exp, coord, *, fit_position=False, neighbours=(), systematic_fraction=0.03):
    """Fit flux/background; free position only when initialized from an independent detection.

    neighbours is an explicit sequence of native (x,y) positions fitted simultaneously.
    No neighbour position, masked pixel value or catalogue flux is invented.
    """
    if not np.isfinite(systematic_fraction) or systematic_fraction < 0:
        raise ValueError("systematic_fraction must be finite and nonnegative")
    x0, y0 = exp.sky_to_pixel(coord)
    if not exp.contains(coord):
        return None
    wl, bw = exp.wavelength_at(x0, y0)
    point = {"wave_um": float(wl) if np.isfinite(wl) else None,
             "bw_um": float(bw) if np.isfinite(bw) else None,
             "flux_ujy": None, "err_ujy": None, "stat_err_ujy": None,
             "mjd": exp.mjd, "det": exp.detector, "obs_id": exp.obs_id,
             "wavelength_source": exp.wavelength_source, "time_source": exp.time_source,
             "time_assumptions": exp.time_assumptions, "provenance": exp.provenance,
             "method": "native_epsf_fit", "quality_reasons": [], "flags": [],
             "position_assisted": not fit_position}
    try:
        model, ovs, zone = select_epsf(exp, x0, y0)
    except ValueError as error:
        point["quality_reasons"].append(str(error))
        return photometry_quality(exp, point, psf=True)
    iy, ix = int(round(y0)), int(round(x0))
    ya, yb = max(0, iy-5), min(exp.shape[0], iy+6)
    xa, xb = max(0, ix-5), min(exp.shape[1], ix+6)
    yy, xx = np.mgrid[ya:yb, xa:xb]
    support = np.hypot(xx-x0, yy-y0) <= 2.5
    region = np.s_[ya:yb, xa:xb]
    mask = exp.bad_mask(include_contamination=True)[region]
    point["flags"] = sorted({name for y,x in zip(*np.nonzero(support))
                             for name in exp.flag_names_at(x+xa, y+ya)})
    if not exp.contains(coord, margin_px=5):
        point["quality_reasons"].append("psf_region_outside_cutout")
    # Saturated cores and contamination are conservatively rejected; fitting
    # remaining pixels still produces a diagnostic result, never a core fill.
    if np.any(mask & support):
        point["quality_reasons"].append("masked_psf_core")
    known = sum(1 << bit for bit in exp.flag_bits.values())
    unknown = (exp.flags[region].astype(np.int64) & ~known) != 0
    if np.any(unknown):
        point["quality_reasons"].append("unknown_pixel_flags")
        mask |= unknown
    # SOURCE is retained throughout the target/neighbor PSFs, and excluded only
    # in background pixels outside their support.
    background = np.hypot(xx-x0,yy-y0)>5
    for nx,ny in neighbours:
        background &= np.hypot(xx-nx,yy-ny)>5
    mask |= background & exp.flag_mask(("SOURCE",))[region]
    # Fit in surface-brightness units. F*PSF/(1e12*omega) has units MJy/sr.
    if exp.calibration is None or exp.calibration.solid_angle is None:
        point["quality_reasons"].append("missing_solid_angle")
        return photometry_quality(exp, point, psf=True)
    omega = exp.calibration.omega_sr(exp, xx, yy)
    mask |= ~np.isfinite(omega) | (omega <= 0)
    good = ~mask
    if good.sum() < 12 + len(neighbours):
        point["quality_reasons"].append("insufficient_psf_pixels")
        return photometry_quality(exp, point, psf=True)
    image, sigma = exp.image[region].astype(float), np.sqrt(exp.variance[region])
    templates = []
    for nx, ny in neighbours:
        nm, no, _ = select_epsf(exp, nx, ny)
        templates.append(render_epsf(nm, no, xx, yy, nx, ny) / (1e12 * omega))
    def basis(x,y):
        return np.column_stack([a[good] for a in
            [render_epsf(model, ovs, xx, yy, x, y) / (1e12*omega), *templates,
             np.ones(xx.shape), xx-x0, yy-y0]])
    A = basis(x0,y0)
    Aw, data = A / sigma[good,None], image[good] / sigma[good]
    initial, *_ = np.linalg.lstsq(Aw, data, rcond=None)
    nflux = 1 + len(neighbours)
    def residual(p):
        b = basis(p[-2],p[-1]) if fit_position else A
        return (b @ p[:A.shape[1]] - image[good]) / sigma[good]
    start = np.r_[initial, x0,y0] if fit_position else initial
    lower, upper = np.full(len(start), -np.inf), np.full(len(start), np.inf)
    lower[:nflux] = 0
    if fit_position:
        lower[-2:], upper[-2:] = [x0-1,y0-1], [x0+1,y0+1]
    start[:nflux] = np.maximum(start[:nflux],1e-9)
    fit = least_squares(residual, start, bounds=(lower,upper), x_scale="jac")
    dof = good.sum()-len(start)
    column_scale = np.linalg.norm(fit.jac,axis=0)
    J = fit.jac / np.maximum(column_scale,1e-30)
    rank = np.linalg.matrix_rank(J)
    cov = np.linalg.pinv(J.T @ J) / np.maximum(np.outer(column_scale,column_scale),1e-60)
    stat = float(np.sqrt(max(cov[0,0],0)))
    flux = float(fit.x[0])
    sys = abs(flux)*systematic_fraction
    x,y = fit.x[-2:] if fit_position else (x0,y0)
    if fit_position and select_epsf(exp,x,y)[2] != zone:
        point["quality_reasons"].append("epsf_zone_boundary_crossed")
    r = residual(fit.x)
    chisq = float(r @ r)
    ql = float(np.mean(np.abs(r[support[good]]))) if support[good].any() else None
    fit_valid = bool(fit.success and dof > 0 and rank == len(start) and chi2.sf(chisq,dof) >= .001
                     and ql is not None and ql <= 2 and not np.any(fit.active_mask))
    # A field detection inside the fit window that is not modelled is a blend.
    from .centroid import detect_all
    detections = detect_all(np.where(exp.bad_mask(True),np.nan, exp.background_subtracted()))
    extra = [c for c in detections
             if 1.5 < np.hypot(c["xy"][0]-x0,c["xy"][1]-y0) < 5
             and not any(np.hypot(c["xy"][0]-nx,c["xy"][1]-ny)<1.5 for nx,ny in neighbours)]
    if extra:
        point["quality_reasons"].append("unmodelled_neighbour")
    point.update(flux_ujy=flux, stat_err_ujy=stat, systematic_err_ujy=sys,
                 err_ujy=float(np.hypot(stat,sys)), systematic_fraction=systematic_fraction,
                 systematic_source="project_assumption", epsf_zone=zone,
                 xy=[float(x),float(y)], covariance_xy_px2=cov[-2:,-2:].tolist() if fit_position else None,
                 fit_chi2=chisq, fit_dof=int(dof), fit_ql=ql, fit_p_value=float(chi2.sf(chisq,dof)) if dof>0 else None,
                 neighbour_count=len(neighbours), n_fit_pixels=int(good.sum()),
                 calibration_provenance=exp.calibration.provenance)
    return photometry_quality(exp, point, psf=True, fit_valid=fit_valid)
