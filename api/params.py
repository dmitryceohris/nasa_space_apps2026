"""Parse and validate query parameters (names and ranges from the API spec, section 2)."""
import warnings
from dataclasses import dataclass

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates.angles import IllegalMinuteWarning, IllegalSecondWarning

from .errors import ApiError

MODES = ("single", "coadd", "diff", "trail")


@dataclass(frozen=True)
class Params:
    target: str | None
    ra: float | None  # degrees, set only when ra+dec were given
    dec: float | None
    fov: float  # arcmin
    det: int
    mode: str
    max_per_epoch: int


def _bad(msg):
    return ApiError(400, "bad_param", msg)


def _num(name, val, lo, hi, cast, default):
    if val in (None, ""):
        return default
    try:
        v = cast(val)
    except (TypeError, ValueError):
        raise _bad(f"{name} должен быть числом")
    if not lo <= v <= hi:
        raise _bad(f"{name} должен быть в диапазоне {lo}..{hi}")
    return v


def parse_coords(ra: str, dec: str):
    """Degrees, or sexagesimal ('08h55m10.83s' / '08:55:10.83')."""
    ra, dec = str(ra).strip(), str(dec).strip()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", IllegalMinuteWarning)
            warnings.simplefilter("error", IllegalSecondWarning)
            warnings.filterwarnings("ignore", message="'hour' was found.*")
            try:
                ra_d = float(ra)
            except ValueError:
                ra_d = float(Angle(ra, unit=u.hourangle).degree)
            try:
                dec_d = float(dec)
            except ValueError:
                dec_d = float(Angle(dec, unit=u.deg).degree)
    except Exception:
        raise _bad("не удалось разобрать ra/dec")
    if not 0 <= ra_d <= 360:
        raise _bad("ra должен быть в диапазоне 0..360")
    if not -90 <= dec_d <= 90:
        raise _bad("dec должен быть в диапазоне -90..90")
    return (0.0 if abs(ra_d - 360) < 1e-12 else ra_d), dec_d


def parse(target=None, ra=None, dec=None, fov=None, det=None, mode=None, max_per_epoch=None) -> Params:
    has_ra, has_dec = ra not in (None, ""), dec not in (None, "")
    if has_ra != has_dec:
        raise _bad("ra и dec нужно передавать вместе")
    if not has_ra and not (target or "").strip():
        raise _bad("нужен либо target, либо пара ra+dec")
    ra_d = dec_d = None
    if has_ra:  # ra/dec win over target
        ra_d, dec_d = parse_coords(ra, dec)
    mode = (mode or "coadd").strip()
    if mode not in MODES:
        raise _bad(f"mode должен быть одним из: {', '.join(MODES)}")
    return Params(
        target=(target or "").strip() or None, ra=ra_d, dec=dec_d,
        fov=_num("fov", fov, 1, 20, float, 5.0),
        det=_num("det", det, 1, 6, int, 3),
        mode=mode,
        max_per_epoch=_num("max_per_epoch", max_per_epoch, 1, 20, int, 6),
    )
