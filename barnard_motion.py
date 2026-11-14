"""Offline barycentric ICRS propagation of Gaia DR3 Barnard's Star.

No annual parallax or observer correction is applied. See BARNARD_MOTION.md.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
import re
from typing import Iterable

import astropy.units as u
from astropy.coordinates import Distance, SkyCoord
from astropy.time import Time
from astropy.utils import iers
import numpy as np

# UTC <-> TCB uses bundled ERFA/leap-second data. No network access is needed.
iers.conf.auto_download = False
iers.conf.auto_max_age = None

SOURCE_ID = "4472832130942575872"
DEFAULT_CATALOG = Path(__file__).resolve().with_name("OBJECT_PARAMETERS.csv")
PARAMETERS = ("ra", "dec", "parallax", "pmra", "pmdec")
VALUE_COLUMNS = ("ra_deg", "dec_deg", "parallax_mas", "pmra_mas_yr",
                 "pmdec_mas_yr", "radial_velocity_km_s")
ERROR_COLUMNS = ("ra_error_mas", "dec_error_mas", "parallax_error_mas",
                 "pmra_error_mas_yr", "pmdec_error_mas_yr",
                 "radial_velocity_error_km_s")
CORRELATIONS = tuple(f"{a}_{b}_corr" for a, b in combinations(PARAMETERS, 2))
TIME_FORMATS = ("iso", "isot", "mjd", "jyear")
TIME_SCALES = ("utc", "tai", "tt", "tcg", "tdb", "tcb")
MAS_PER_DEG = 3_600_000.0
# Central differences: input units are mas, mas, mas, mas/yr, mas/yr, km/s.
# Steps exceed floating-point angular noise, while remaining locally linear.
JACOBIAN_STEPS = np.array([10.0, 10.0, 0.1, 0.1, 0.1, 0.1])
OUTPUT_COLUMNS = (
    "obs_id", "detector", "data_release", "access_url", "is_demo",
    "input_time_format", "input_time_scale", "input_time", "input_t_min",
    "input_t_max", "time_method", "time_utc", "time_tcb", "time_mjd_utc",
    "time_jyear_tcb", "coordinate_frame", "position_type", "ra_deg", "dec_deg",
    "sigma_alpha_star_mas", "sigma_delta_mas", "ra_dec_corr",
    "cov_alpha_star_delta_mas2", "distance_pc", "parallax_scale_1au_mas",
    "parallax_applied", "source_id", "gaia_source_url",
)


@dataclass(frozen=True)
class GaiaSource:
    epoch: Time
    values: np.ndarray  # RA/Dec in deg; parallax in mas; PM in mas/yr; RV in km/s.
    covariance: np.ndarray  # First two axes are alpha-star and delta, in mas.
    source_url: str
    row: dict[str, str]
    source_id: str = SOURCE_ID


@dataclass(frozen=True)
class Observation:
    time: Time
    obs_id: str = ""
    metadata: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class Prediction:
    coordinate: SkyCoord
    covariance_mas2: np.ndarray  # 2x2 local east/north covariance at output epoch.
    jacobian: np.ndarray  # 2x6; output angular axes in mas.


def _read_csv(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        names = reader.fieldnames or []
        if not names or len(names) != len(set(names)):
            raise ValueError(f"{path}: отсутствует заголовок или повторяются имена колонок")
        rows = list(reader)
    for number, row in enumerate(rows, 2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"{path}, строка {number}: число полей не совпадает с заголовком")
    return names, rows


def _number(row: dict[str, str], key: str) -> float:
    try:
        result = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Параметр {key}: требуется заполненное числовое значение") from exc
    if not np.isfinite(result):
        raise ValueError(f"Параметр {key}: значение должно быть конечным")
    return result


def load_source(path: str | Path = DEFAULT_CATALOG) -> GaiaSource:
    """Load only the exact DR3 identifier, preserving IDs as strings."""
    names, rows = _read_csv(path)
    required = {"source_id", "canonical_id", "catalog", "coordinate_frame", "ref_epoch",
                "ref_epoch_format", "ref_epoch_scale", "ra_error_basis", "source_url",
                *VALUE_COLUMNS, *ERROR_COLUMNS, *CORRELATIONS}
    missing = sorted(required - set(names))
    if missing:
        raise ValueError("В таблице параметров отсутствуют колонки: " + ", ".join(missing))
    selected = [row for row in rows if row["source_id"] == SOURCE_ID]
    if len(selected) != 1:
        raise ValueError(f"Ожидается ровно одна запись Gaia DR3 {SOURCE_ID}; найдено {len(selected)}")
    row = selected[0]
    if row["canonical_id"] != f"Gaia DR3 {SOURCE_ID}" or "Gaia DR3" not in row["catalog"]:
        raise ValueError("Идентификатор должен относиться именно к Gaia DR3")
    if (row["coordinate_frame"] != "ICRS" or row["ref_epoch_format"] != "jyear"
            or row["ref_epoch_scale"].lower() != "tcb" or _number(row, "ref_epoch") != 2016.0):
        raise ValueError("Для этой записи требуются ICRS и эпоха J2016.0 TCB")
    if row["ra_error_basis"] != "sigma_alpha_cosdec":
        raise ValueError("Ошибка RA должна быть задана как sigma_alpha_cosdec (σα⋆)")
    if not row["source_url"].startswith("https://"):
        raise ValueError("Параметр source_url: требуется ссылка на исходную запись Gaia")
    values = np.array([_number(row, key) for key in VALUE_COLUMNS])
    errors = np.array([_number(row, key) for key in ERROR_COLUMNS])
    if not (0 <= values[0] < 360 and -90 < values[1] < 90 and values[2] > 0):
        raise ValueError("Требуются RA в [0, 360), Dec в (-90, 90) и положительный параллакс")
    if np.any(errors <= 0):
        raise ValueError("Все формальные погрешности, включая RV, должны быть положительными")
    correlation = np.eye(6)
    for i, j in combinations(range(5), 2):
        key = f"{PARAMETERS[i]}_{PARAMETERS[j]}_corr"
        rho = _number(row, key)
        if not -1 <= rho <= 1:
            raise ValueError(f"Параметр {key}: корреляция должна лежать в [-1, 1]")
        correlation[i, j] = correlation[j, i] = rho
    if np.linalg.eigvalsh(correlation).min() < -1e-12:
        raise ValueError("Матрица корреляций Gaia не является положительно полуопределённой")
    covariance = correlation * np.outer(errors, errors)
    return GaiaSource(Time(2016.0, format="jyear", scale="tcb"), values,
                      covariance, row["source_url"], row)


def parse_time(value: str | float, time_format: str | None = None,
               time_scale: str | None = None) -> Time:
    """Parse one date. Numeric dates require explicit format AND scale."""
    fmt = (time_format or "iso").lower()
    if fmt not in TIME_FORMATS:
        raise ValueError(f"Неизвестный формат времени {fmt!r}; доступны {', '.join(TIME_FORMATS)}")
    if fmt in ("mjd", "jyear") and time_scale is None:
        raise ValueError("Для числовой даты укажите формат (mjd/jyear) и шкалу времени явно")
    scale = (time_scale or "utc").lower()
    if scale not in TIME_SCALES:
        raise ValueError(f"Неизвестная или неподдерживаемая шкала {scale!r}; доступны "
                         f"{', '.join(TIME_SCALES)}. UT1 требует отдельного учёта вращения Земли")
    text = str(value).strip()
    if fmt in ("iso", "isot"):
        pattern = (r"\d{4}-\d{2}-\d{2}(?:[T ](?:[01]\d|2[0-3]):[0-5]\d:"
                   r"(?:[0-5]\d|60)(?:\.\d+)?)?Z?")
        if not re.fullmatch(pattern, text):
            raise ValueError("ISO-дата: ожидается YYYY-MM-DD или YYYY-MM-DDTHH:MM:SS; "
                             "числовые даты требуют явных --format и --scale; смещения часового пояса не принимаются")
        if text.endswith("Z"):
            if scale != "utc":
                raise ValueError("Суффикс Z допустим только для шкалы UTC")
            text = text[:-1]
        text = text.replace(" ", "T")
        if len(text) == 10:
            text += "T00:00:00"
        fmt = "isot"
    else:
        try:
            if not np.isfinite(float(text)):
                raise ValueError
        except ValueError as exc:
            raise ValueError(f"Дата в формате {fmt}: требуется конечное число") from exc
    try:
        result = Time(text, format=fmt, scale=scale, precision=9)
        # Reject normalised/nonexistent leap seconds rather than silently shifting dates.
        if fmt == "isot" and text[17:19] == "60" and result.isot[17:19] != "60":
            raise ValueError("на эту дату/шкалу дополнительная секунда не существует")
        if not result.isscalar or not np.isfinite(result.jd):
            raise ValueError("требуется один конечный момент времени")
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"Некорректная дата {value!r} ({fmt}, {scale}): {exc}") from exc
    return result


def observation_from_row(row: dict[str, str]) -> Observation:
    """Validate an exact time or interval; midpoint is in continuous TAI."""
    fmt, scale = row.get("time_format", "").strip(), row.get("time_scale", "").strip()
    if not fmt or not scale:
        raise ValueError("Обязательны time_format и time_scale в каждой строке")
    exact, start, end = (row.get(key, "").strip() for key in ("time", "t_min", "t_max"))
    if exact and (start or end):
        raise ValueError("Неоднозначная строка: задайте либо time, либо пару t_min/t_max")
    if exact:
        time, method = parse_time(exact, fmt, scale), "exact"
    elif start and end:
        t0, t1 = parse_time(start, fmt, scale), parse_time(end, fmt, scale)
        if t1 < t0:
            raise ValueError("Обратный интервал: t_max раньше t_min")
        time, method = t0.tai + (t1.tai - t0.tai) / 2, "midpoint_tai"
    else:
        raise ValueError("Требуется time или полная пара t_min и t_max")
    demo = row.get("is_demo", "false").strip().lower() or "false"
    if demo not in ("true", "false"):
        raise ValueError("Поле is_demo должно быть true или false")
    obs_id = row.get("obs_id", "")
    if not obs_id.strip() and demo != "true":
        raise ValueError("Для реального наблюдения нужен obs_id; учебная строка требует is_demo=true")
    metadata = {key: row.get(key, "") for key in ("detector", "data_release", "access_url")}
    metadata.update(is_demo=demo, input_time_format=fmt.lower(), input_time_scale=scale.lower(),
                    input_time=exact, input_t_min=start, input_t_max=end, time_method=method)
    return Observation(time, obs_id, metadata)


def load_observations(path: str | Path) -> list[Observation]:
    names, rows = _read_csv(path)
    missing = {"obs_id", "time_format", "time_scale"} - set(names)
    if missing:
        raise ValueError("В CSV наблюдений отсутствуют колонки: " + ", ".join(sorted(missing)))
    if not rows:
        raise ValueError("CSV наблюдений не содержит строк")
    observations = []
    for number, row in enumerate(rows, 2):
        try:
            observations.append(observation_from_row(row))
        except ValueError as exc:
            raise ValueError(f"{path}, строка {number}, obs_id={row.get('obs_id', '')!r}: {exc}") from exc
    return observations


def coordinate_at(source: GaiaSource, time: Time,
                  perturbation: np.ndarray | None = None) -> SkyCoord:
    """Evaluate the requested Astropy space-motion model, optionally perturbed."""
    if not isinstance(time, Time) or not time.isscalar:
        raise ValueError("Передайте один момент astropy.time.Time")
    values = source.values.copy()
    if perturbation is not None:
        delta = np.asarray(perturbation, dtype=float)
        if delta.shape != (6,) or not np.isfinite(delta).all():
            raise ValueError("Возмущение должно содержать шесть конечных компонент")
        values[0] += delta[0] / (MAS_PER_DEG * np.cos(np.deg2rad(source.values[1])))
        values[1] += delta[1] / MAS_PER_DEG
        values[2:] += delta[2:]
    ra, dec, plx, pmra, pmdec, rv = values
    initial = SkyCoord(ra=ra * u.deg, dec=dec * u.deg,
                       distance=Distance(parallax=plx * u.mas),
                       pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                       radial_velocity=rv * u.km / u.s, obstime=source.epoch, frame="icrs")
    return initial.apply_space_motion(new_obstime=time.tcb)


def tangent_offsets_mas(reference: SkyCoord, other: SkyCoord) -> np.ndarray:
    """Local east/north spherical offsets; differential basis is dRA*cos(Dec), dDec."""
    # Position-only coordinates avoid comparing frames with different obstime metadata.
    origin = SkyCoord(ra=reference.ra, dec=reference.dec, frame="icrs")
    target = SkyCoord(ra=other.ra, dec=other.dec, frame="icrs")
    east, north = origin.spherical_offsets_to(target)
    return np.array([east.to_value(u.mas), north.to_value(u.mas)])


def predict(source: GaiaSource, time: Time, *, step_scale: float = 1.0) -> Prediction:
    """Propagate the formal 6D covariance with a numerical 2x6 Jacobian."""
    if not np.isfinite(step_scale) or step_scale <= 0:
        raise ValueError("Множитель шага Якобиана должен быть положительным")
    center = coordinate_at(source, time)
    jacobian = np.empty((2, 6))
    for axis, step in enumerate(JACOBIAN_STEPS * step_scale):
        delta = np.zeros(6)
        delta[axis] = step
        plus = tangent_offsets_mas(center, coordinate_at(source, time, delta))
        minus = tangent_offsets_mas(center, coordinate_at(source, time, -delta))
        jacobian[:, axis] = (plus - minus) / (2 * step)
    covariance = jacobian @ source.covariance @ jacobian.T
    covariance = (covariance + covariance.T) / 2
    return Prediction(center, covariance, jacobian)


def calculate_observation(source: GaiaSource, observation: Observation) -> dict[str, object]:
    prediction = predict(source, observation.time)
    coord, covariance = prediction.coordinate, prediction.covariance_mas2
    sigma = np.sqrt(np.diag(covariance))
    time = observation.time.replicate()
    time.precision = 9
    row: dict[str, object] = {key: "" for key in OUTPUT_COLUMNS}
    row.update(observation.metadata)
    row.update(obs_id=observation.obs_id, time_utc=time.utc.isot + "Z", time_tcb=time.tcb.isot,
               time_mjd_utc=f"{time.utc.mjd:.12f}", time_jyear_tcb=f"{time.tcb.jyear:.12f}",
               coordinate_frame="ICRS", position_type="barycentric_astrometric",
               ra_deg=float(coord.ra.deg), dec_deg=float(coord.dec.deg),
               sigma_alpha_star_mas=float(sigma[0]), sigma_delta_mas=float(sigma[1]),
               ra_dec_corr=float(np.clip(covariance[0, 1] / np.prod(sigma), -1, 1)),
               cov_alpha_star_delta_mas2=float(covariance[0, 1]),
               distance_pc=float(coord.distance.to_value(u.pc)),
               parallax_scale_1au_mas=float(coord.distance.to(u.mas, equivalencies=u.parallax()).value),
               parallax_applied="false", source_id=source.source_id, gaia_source_url=source.source_url)
    return row


def calculate_observations(source: GaiaSource, observations: Iterable[Observation]) -> list[dict[str, object]]:
    return [calculate_observation(source, observation) for observation in observations]


def write_results(rows: Iterable[dict[str, object]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Барицентрические астрометрические координаты Барнарда ICRS; без годичного параллакса")
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--date", help="ISO (по умолчанию UTC), MJD или юлианский год")
    mode.add_argument("--observations", type=Path, help="CSV с time_format/time_scale и time либо t_min/t_max")
    parser.add_argument("--format", choices=TIME_FORMATS, dest="time_format")
    parser.add_argument("--scale", dest="time_scale", help="utc, tai, tt, tcg, tdb, tcb")
    parser.add_argument("--output", type=Path, default=Path("output/barnard_positions.csv"))
    args = parser.parse_args(argv)
    try:
        source = load_source(args.catalog)
        if args.observations:
            if args.time_format or args.time_scale:
                raise ValueError("Для --observations формат и шкала задаются в CSV, без --format/--scale")
            observations = load_observations(args.observations)
        else:
            time = parse_time(args.date, args.time_format, args.time_scale)
            observations = [Observation(time, metadata={"input_time": args.date,
                            "input_time_format": args.time_format or "iso",
                            "input_time_scale": args.time_scale or "utc", "time_method": "exact",
                            "is_demo": "false"})]
        rows = calculate_observations(source, observations)
        write_results(rows, args.output)
    except (ValueError, OSError, csv.Error) as exc:
        parser.exit(2, f"Ошибка: {exc}\n")
    print(f"Сохранено {len(rows)} строк: {args.output}")
    print("ICRS, барицентрические астрометрические координаты; годичный параллакс не применён.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
