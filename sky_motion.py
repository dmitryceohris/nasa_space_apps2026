"""Five-target space motion and geometric annual parallax; see docs/PHYSICS.md.

Only --prepare-observations accesses the network. No photometry or detection is
performed. The original Barnard implementation remains the space-motion engine.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
from io import BytesIO
import json
from pathlib import Path
import ssl
import time as wall_time
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import astropy.units as u
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.io import fits
from astropy.io.votable import parse
from astropy.time import Time
from astropy.wcs import WCS
import numpy as np

from barnard_motion import DEFAULT_CATALOG, VALUE_COLUMNS, coordinate_at, parse_time

ROOT = Path(__file__).resolve().parent
TARGET_IDS = {
    "barnard": "Barnard's Star", "epsilon_indi_b": "Epsilon Indi B",
    "wise_0855": "WISE 0855-0714", "rw_aur": "RW Aur", "r_dra": "R Dra",
}
PREDICTION_COLUMNS = (
    "target", "epoch_id", "obs_id", "detector", "release", "mjd", "time_source",
    "ra_bary", "dec_bary", "ra_geo", "dec_geo", "d_east_arcsec",
    "d_north_arcsec", "shift_px", "coverage_margin_arcsec", "quality_status", "assumptions",
)
OBSERVATION_COLUMNS = (
    "target", "epoch_id", "obs_id", "detector", "release", "mjd", "time_source",
    "coverage_margin_arcsec", "quality_status", "assumptions", "access_url", "cutout_url",
)
COMMON_ASSUMPTIONS = ("geocenter_geometric", "no_aberration", "no_new_covariance")
BAD_FLAGS = {"OVERFLOW", "NONLINEAR", "SUR_ERROR", "NONFUNC", "MISSING_DATA", "HOT", "COLD"}


@dataclass(frozen=True)
class Target:
    name: str
    epoch: Time
    values: np.ndarray
    row: dict[str, str]
    assumptions: tuple[str, ...]


@dataclass(frozen=True)
class Position:
    barycentric: SkyCoord
    geocentric: SkyCoord
    assumptions: tuple[str, ...]


def load_targets(path: str | Path = DEFAULT_CATALOG) -> dict[str, Target]:
    """Read all five solutions without inventing missing measurements/errors."""
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    targets = {}
    for row in rows:
        name = row["object"]
        if name in targets or row["coordinate_frame"] != "ICRS":
            raise ValueError(f"Duplicate target or unsupported frame: {name}")
        assumptions = list(COMMON_ASSUMPTIONS)
        scale = row["ref_epoch_scale"].strip().lower()
        if not scale:
            if name != TARGET_IDS["wise_0855"]:
                raise ValueError(f"Unknown reference epoch scale: {name}")
            scale = "utc"
            assumptions.append("reference_epoch_scale_assumed_utc")
        epoch = parse_time(row["ref_epoch"], row["ref_epoch_format"], scale)
        values = [float(row[key]) for key in VALUE_COLUMNS[:-1]]
        rv = row["radial_velocity_km_s"].strip()
        values.append(float(rv) if rv else 0.0)
        if not rv:
            assumptions.append("rv_assumed_zero")
        if not np.isfinite(values).all() or not (0 <= values[0] < 360 and -90 < values[1] < 90 and values[2] > 0):
            raise ValueError(f"Invalid astrometric parameters: {name}")
        if row.get("ruwe") and float(row["ruwe"]) > 1.4:
            assumptions.append("gaia_ruwe_high")
        if name == TARGET_IDS["epsilon_indi_b"]:
            assumptions.append("Ba_position_not_BaBb_photocentre")
        if name == TARGET_IDS["rw_aur"]:
            assumptions.append("literature_A_eastern_component_blended_pair")
        if name == TARGET_IDS["wise_0855"]:
            assumptions.append("position_error_basis_and_correlations_unpublished")
        targets[name] = Target(name, epoch, np.array(values), row, tuple(assumptions))
    if set(targets) != set(TARGET_IDS.values()):
        raise ValueError("Catalogue must contain exactly the five challenge targets")
    return targets


def predict_position(target: Target | str, time: Time) -> Position:
    """Barycentric and observer-geometric directions in ICRS-aligned axes.

    The geocentric result is NOT a GCRS apparent coordinate: no aberration,
    light deflection, spacecraft orbit or additional light-time correction.
    """
    if isinstance(target, str):
        target = load_targets()[TARGET_IDS.get(target, target)]
    barycentric = coordinate_at(target, time)
    earth = get_body_barycentric("earth", time, ephemeris="builtin")
    geocentric = SkyCoord(barycentric.cartesian.without_differentials() - earth, frame="icrs")
    return Position(barycentric, geocentric, target.assumptions)


def angular_position(coord: SkyCoord) -> SkyCoord:
    return SkyCoord(ra=coord.ra, dec=coord.dec, frame="icrs")


def footprint_margin(region: str, coord: SkyCoord) -> float:
    """Signed distance to edges of a convex spherical SIA POLYGON, arcsec."""
    parts = region.split()
    if len(parts) < 8 or parts[:2] != ["POLYGON", "ICRS"]:
        raise ValueError("Expected a convex POLYGON ICRS footprint")
    vertices = np.array([float(value) for value in parts[2:]]).reshape(-1, 2)
    vectors = SkyCoord(vertices[:, 0] * u.deg, vertices[:, 1] * u.deg).cartesian.xyz.value.T
    if np.allclose(vectors[0], vectors[-1], atol=1e-12, rtol=0):
        vectors = vectors[:-1]
    normals = np.cross(vectors, np.roll(vectors, -1, axis=0))
    norms = np.linalg.norm(normals, axis=1)
    if len(vectors) < 3 or np.any(norms < 1e-14):
        raise ValueError("Degenerate footprint")
    normals /= norms[:, None]
    normals *= np.where(normals @ vectors.mean(axis=0) >= 0, 1, -1)[:, None]
    distances = np.arcsin(np.clip(normals @ angular_position(coord).cartesian.xyz.value, -1, 1))
    return float(np.rad2deg(distances.min()) * 3600)


def select_epochs(rows: list[dict]) -> list[tuple[str, list[dict]]]:
    """One selection per >30-day-gap cluster, or three continuous windows."""
    ordered = sorted(rows, key=lambda row: float(row["mjd"]))
    if not ordered:
        return []
    groups = [[ordered[0]]]
    for row in ordered[1:]:
        if float(row["mjd"]) - float(groups[-1][-1]["mjd"]) > 30:
            groups.append([])
        groups[-1].append(row)
    prefix = "cluster"
    if len(groups) == 1 and float(ordered[-1]["mjd"]) > float(ordered[0]["mjd"]):
        first, last = float(ordered[0]["mjd"]), float(ordered[-1]["mjd"])
        bins = [[], [], []]
        for row in ordered:
            index = min(2, int(3 * (float(row["mjd"]) - first) / (last - first)))
            bins[index].append(row)
        groups, prefix = [group for group in bins if group], "window"
    if len(groups) > 3:
        groups = [groups[0], groups[len(groups) // 2], groups[-1]]
    return [(f"{prefix}_{i + 1}", group) for i, group in enumerate(groups)]


def observation_time(header: fits.Header, row: dict) -> tuple[Time, str, list[str]]:
    scale = str(header.get("TIMESYS", "")).strip().lower()
    assumptions = [] if scale else ["observation_time_scale_assumed_utc"]
    scale = scale or "utc"
    if "MJD-AVG" in header:
        return parse_time(header["MJD-AVG"], "mjd", scale), "fits_mjd_avg", assumptions
    if "DATE-AVG" in header:
        return parse_time(header["DATE-AVG"], "isot", scale), "fits_date_avg", assumptions
    start = parse_time(row["t_min"], "mjd", "utc")
    end = parse_time(row["t_max"], "mjd", "utc")
    if end < start:
        raise ValueError("Reversed SIA time interval")
    # SIA metadata do not supply TIMESYS; do not infer its scale from an unrelated header.
    return start.tai + (end.tai - start.tai) / 2, "sia_midpoint", ["sia_time_scale_assumed_utc"]


def assess_flags(data: np.ndarray, header: fits.Header, x: float, y: float) -> tuple[str, dict]:
    """Inspect the local 5x5 region using definitions supplied by this product."""
    if data.ndim != 2 or not np.isfinite([x, y]).all():
        raise ValueError("Expected a 2D FLAGS image and finite target pixel")
    ix, iy = int(np.floor(x + 0.5)), int(np.floor(y + 0.5))
    if ix < 2 or iy < 2 or ix + 2 >= data.shape[1] or iy + 2 >= data.shape[0]:
        return "unchecked", {"reason": "incomplete_5x5_flags_region"}
    definitions = {int(key[4:]): str(header[key]).strip() for key in header if key.startswith("MSKN") and key[4:].isdigit()}
    known = set(definitions.values())
    if not BAD_FLAGS.issubset(known):
        return "unchecked", {"reason": "missing_flag_definitions"}
    patch = np.asarray(data[iy - 2:iy + 3, ix - 2:ix + 3], dtype=np.uint64)
    hits = {name: int(np.count_nonzero(patch & (np.uint64(1) << np.uint64(bit))))
            for bit, name in definitions.items() if name in BAD_FLAGS}
    evidence = {"pixel_x": x, "pixel_y": y, "flags_5x5": patch.tolist(), "bad_flag_counts": hits,
                "definitions": definitions}
    return ("flagged" if any(hits.values()) else "usable"), evidence


def fetch(url: str, *, max_bytes: int = 8_000_000) -> bytes:
    """Verified TLS, 15-second socket timeout, one retry, bounded response size."""
    try:
        import certifi  # Optional existing environment CA bundle; never disable TLS verification.
    except ImportError:
        context = ssl.create_default_context()
    else:
        context = ssl.create_default_context(cafile=certifi.where())
    for attempt in range(2):
        try:
            with urlopen(url, timeout=15, context=context) as response:
                if int(response.headers.get("Content-Length", 0)) > max_bytes:
                    raise ValueError("Response exceeds size limit; full images are not downloaded")
                result = response.read(max_bytes + 1)
            if len(result) > max_bytes:
                raise ValueError("Response exceeds size limit")
            return result
        except (URLError, OSError, TimeoutError):
            if attempt:
                raise
    raise RuntimeError("Unreachable")


def query_observations(target: Target, directory: Path) -> tuple[list[dict], dict]:
    query = urlencode([
        ("COLLECTION", "spherex_dr1"), ("COLLECTION", "spherex_qr3"),
        ("POS", f"CIRCLE {target.values[0]} {target.values[1]} 0.05"),
        ("DPTYPE", "image"), ("CALIB", "2"), ("MAXREC", "10000"), ("RESPONSEFORMAT", "VOTABLE"),
    ])
    url = "https://irsa.ipac.caltech.edu/SIA?" + query
    data = fetch(url)
    (directory / "search_results.vot").write_bytes(data)
    votable = parse(BytesIO(data))
    resources = [resource for resource in votable.resources if resource.type == "results"]
    statuses = [str(info.value) for resource in resources for info in resource.infos if info.name == "QUERY_STATUS"]
    if not statuses or any(status != "OK" for status in statuses):
        raise ValueError(f"Incomplete/failed SIA response: {statuses}")
    table = next(table for resource in resources for table in resource.tables).to_table(use_names_over_ids=True)
    rows = []
    for record in table:
        row = {key: "" if np.ma.is_masked(record[key]) else str(record[key]) for key in table.colnames}
        moment, method, assumptions = observation_time(fits.Header(), row)
        row.update(mjd=float(moment.utc.mjd), time_source=method, assumptions=";".join(assumptions))
        rows.append(row)
    return rows, {"query_url": url, "candidates": len(rows), "query_status": statuses}


def check_candidate(target: Target, candidate: dict, directory: Path, epoch_id: str, attempt: int) -> tuple[dict, dict]:
    row = dict(candidate)
    moment, method, assumptions = observation_time(fits.Header(), row)
    position = predict_position(target, moment)
    margin = footprint_margin(row["s_region"], position.geocentric)
    result = {"target": target.name, "epoch_id": epoch_id, "obs_id": row["obs_id"],
              "detector": row["energy_bandpassname"], "release": row["obs_collection"],
              "mjd": float(moment.utc.mjd), "time_source": method, "coverage_margin_arcsec": margin,
              "quality_status": "no_coverage" if margin <= 10 else "unchecked",
              "assumptions": ";".join(assumptions), "access_url": row["access_url"], "cutout_url": ""}
    evidence = {"candidate": row, "result": result}
    if margin <= 10:
        return result, evidence
    access_url = row["access_url"]
    if "image/fits" not in row.get("access_format", "").lower():
        result["assumptions"] += ";unsupported_access_format"
        return result, evidence
    cutout_url = access_url + ("&" if "?" in access_url else "?") + urlencode({
        "center": f"{position.geocentric.ra.deg:.12f},{position.geocentric.dec.deg:.12f}", "size": "0.025"})
    result["cutout_url"] = cutout_url
    try:
        data = fetch(cutout_url)
        with fits.open(BytesIO(data), memmap=False) as hdus:
            image = hdus["IMAGE"]
            if max(image.data.shape) > 128:
                raise ValueError("Server did not return a small cutout")
            stem = f"{epoch_id}_attempt{attempt}"
            for name in ("IMAGE", "FLAGS"):
                hdus[name].header.tofile(directory / f"{stem}_{name.lower()}_header.txt", sep="\n", endcard=True, padding=False, overwrite=True)
            moment, method, assumptions = observation_time(image.header, row)
            position = predict_position(target, moment)
            margin = footprint_margin(row["s_region"], position.geocentric)
            x, y = WCS(image.header).celestial.world_to_pixel(angular_position(position.geocentric))
            quality, flags = assess_flags(hdus["FLAGS"].data, hdus["FLAGS"].header, float(x), float(y))
            if not (0 <= x < image.data.shape[1] and 0 <= y < image.data.shape[0]) or margin <= 10:
                quality = "no_coverage"
            if not np.isfinite(image.data[max(0, int(round(float(y))) - 2):int(round(float(y))) + 3,
                                          max(0, int(round(float(x))) - 2):int(round(float(x))) + 3]).all() and quality == "usable":
                quality = "flagged"
                flags["nonfinite_image_pixels"] = True
            result.update(mjd=float(moment.utc.mjd), time_source=method, assumptions=";".join(assumptions),
                          coverage_margin_arcsec=margin, quality_status=quality)
            evidence.update(flags=flags, cutout_bytes=len(data), cutout_sha256=hashlib.sha256(data).hexdigest(),
                            cutout_shape=list(image.data.shape), headers_stem=stem)
    except (ValueError, OSError, URLError, KeyError, StopIteration) as exc:
        evidence["error"] = str(exc)
        result["assumptions"] += ";quality_check_failed"
    return result, evidence


def write_csv(rows: list[dict], path: Path, columns: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def choose_attempt(attempts: list[dict]) -> dict:
    """Prefer verified evidence to an unchecked replacement; neither is a detection."""
    priority = {"usable": 0, "flagged": 1, "unchecked": 2, "no_coverage": 3}
    return min((item["result"] for item in attempts), key=lambda row: priority[row["quality_status"]])


def prepare_observations(targets: dict[str, Target], output: Path) -> list[dict]:
    """Bounded archive preparation; save partial progress after every target."""
    deadline = wall_time.monotonic() + 15 * 60
    directory = output.parent / "physics_evidence"
    directory.mkdir(parents=True, exist_ok=True)
    all_rows, report = [], {"checked_utc": Time.now().utc.isot, "targets": {}}
    seeds_path = ROOT / "data/irsa_barnard_20261114/observations.csv"
    with seeds_path.open(newline="", encoding="utf-8") as stream:
        barnard_seeds = {row["obs_id"] for row in csv.DictReader(stream)}
    for alias, name in TARGET_IDS.items():
        target = targets[name]
        target_directory = directory / alias
        target_directory.mkdir(exist_ok=True)
        target_report = {"source_id": target.row["source_id"], "source_url": target.row["source_url"], "epochs": []}
        report["targets"][name] = target_report
        try:
            if wall_time.monotonic() >= deadline:
                raise TimeoutError("Preparation time budget exhausted")
            if alias == "barnard":
                saved = ROOT / "data/irsa_barnard_20261114/search_results.vot"
                votable = parse(saved)
                table = next(table for resource in votable.resources if resource.type == "results" for table in resource.tables).to_table(use_names_over_ids=True)
                rows = [{key: "" if np.ma.is_masked(record[key]) else str(record[key]) for key in table.colnames} for record in table]
                for row in rows:
                    row["mjd"] = float(observation_time(fits.Header(), row)[0].utc.mjd)
                target_report["saved_search"] = str(saved.relative_to(ROOT))
            else:
                rows, query = query_observations(target, target_directory)
                target_report.update(query)
            preferred = "SPHEREx-D6" if alias == "wise_0855" else "SPHEREx-D1"
            for epoch_id, group in select_epochs(rows):
                centre = (min(float(row["mjd"]) for row in group) + max(float(row["mjd"]) for row in group)) / 2
                def rank(row):
                    seed = alias == "barnard" and row["obs_id"] in barnard_seeds and row["energy_bandpassname"] == preferred
                    return (not seed, row["obs_collection"] != "spherex_dr1", row["energy_bandpassname"] != preferred,
                            abs(float(row["mjd"]) - centre), row["obs_id"])
                ranked = sorted(group, key=rank)
                first = ranked[0]
                alternate = next((row for row in ranked[1:] if row["energy_bandpassname"] != first["energy_bandpassname"]), None)
                candidates = [first] + ([alternate] if alternate is not None else ranked[1:2])
                attempts = []
                for i, candidate in enumerate(candidates, 1):
                    if wall_time.monotonic() >= deadline:
                        break
                    result, evidence = check_candidate(target, candidate, target_directory, epoch_id, i)
                    attempts.append(evidence)
                    if result["quality_status"] == "usable":
                        break
                if attempts:
                    chosen = choose_attempt(attempts)
                    all_rows.append(chosen)
                    target_report["epochs"].append({"epoch_id": epoch_id, "attempts": attempts, "selected": chosen})
                    print(f"{name}: {epoch_id}, {chosen['detector']}, {chosen['quality_status']}", flush=True)
        except (ValueError, OSError, URLError, KeyError, StopIteration) as exc:
            target_report["error"] = str(exc)
        if not any(row["target"] == name for row in all_rows):
            all_rows.append({"target": name, "quality_status": "unchecked", "assumptions": "no_observations_prepared"})
        write_csv(all_rows, output, OBSERVATION_COLUMNS)
        (directory / "selection.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return all_rows


def calculate_predictions(observations: list[dict], targets: dict[str, Target]) -> list[dict]:
    rows, first_positions = [], {}
    ordered = sorted(observations, key=lambda row: (row["target"], float(row["mjd"]) if row.get("mjd") else float("inf")))
    for observation in ordered:
        name = observation["target"]
        target = targets[name]
        row = {key: observation.get(key, "") for key in PREDICTION_COLUMNS}
        row["assumptions"] = ";".join(dict.fromkeys([*target.assumptions, *filter(None, row["assumptions"].split(";"))]))
        if observation.get("mjd"):
            moment = parse_time(observation["mjd"], "mjd", "utc")
            position = predict_position(target, moment)
            geo = angular_position(position.geocentric)
            first = first_positions.setdefault(name, geo)
            east, north = first.spherical_offsets_to(geo)
            row.update(ra_bary=position.barycentric.ra.deg, dec_bary=position.barycentric.dec.deg,
                       ra_geo=geo.ra.deg, dec_geo=geo.dec.deg, d_east_arcsec=east.arcsec,
                       d_north_arcsec=north.arcsec, shift_px=first.separation(geo).arcsec / 6.15)
        rows.append(row)
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-observations", action="store_true")
    mode.add_argument("--observations", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        targets = load_targets()
        if args.prepare_observations:
            rows = prepare_observations(targets, args.output)
        else:
            with args.observations.open(newline="", encoding="utf-8-sig") as stream:
                observations = list(csv.DictReader(stream))
            rows = calculate_predictions(observations, targets)
            write_csv(rows, args.output, PREDICTION_COLUMNS)
    except (ValueError, OSError, KeyError, csv.Error) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(f"Saved {len(rows)} rows: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
