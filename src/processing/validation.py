"""One conservative quality contract shared by scientific outputs and diagnostics."""
from __future__ import annotations

import numpy as np


def quality(checks: dict[str, bool | None], rejected: list[str] = ()) -> dict:
    """False/unknown never passes. Rejected names are fatal; unknowns are provisional."""
    reasons = [name for name, passed in checks.items() if passed is not True]
    status = "rejected" if any(checks.get(name) is False for name in rejected) else (
        "provisional" if reasons else "validated")
    return {"status": status, "reasons": reasons, "checks": checks}


def accepted(point: dict) -> bool:
    """Legacy ok=False is always a veto; unclassified legacy points are not validated."""
    return bool(point.get("ok") and point.get("validation", {}).get("status") == "validated"
                and point.get("flux_ujy") is not None and point.get("err_ujy") is not None
                and np.isfinite([point["flux_ujy"], point["err_ujy"]]).all()
                and point["err_ujy"] > 0)


def photometry_quality(exp, point: dict, *, psf: bool = False, fit_valid: bool | None = None) -> dict:
    flux, err = point.get("flux_ujy"), point.get("stat_err_ujy", point.get("err_ujy"))
    finite = flux is not None and err is not None and bool(np.isfinite([flux, err]).all()) and err > 0
    checks = {reason: False for reason in point.get("quality_reasons", [])}
    checks.update(flags_missing=bool(exp.flags_available),
                  unverified_time=bool(exp.time_source in ("MJD-AVG", "DATE-AVG") and not exp.time_assumptions),
                  approximate_wavelength=exp.wavelength_source == "cwave_cband",
                  preliminary_aperture=psf,
                  missing_solid_angle=bool(exp.calibration is not None and exp.calibration.solid_angle is not None),
                  nonfinite_measurement=bool(finite),
                  low_snr=bool(finite and flux / err >= 5),
                  poor_psf_fit=fit_valid if psf else None)
    checks["invalid_wavelength"]=bool(point.get("wave_um") is not None and point.get("bw_um") is not None
        and np.isfinite([point["wave_um"],point["bw_um"]]).all() and point["wave_um"]>0 and point["bw_um"]>0)
    checks["unverified_wcs_frame"]=bool(exp.header is not None and str(exp.header.get("RADESYS","")).strip()=="ICRS")
    fatal = [r for r in point.get("quality_reasons", []) if r != "invalid_wavelength"]
    result = quality(checks, fatal + ["nonfinite_measurement", "low_snr", "poor_psf_fit"])
    point.update(validation=result, ok=result["status"] == "validated", quality_reasons=result["reasons"])
    return point
