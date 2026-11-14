"""Scientific acceptance checks and input/CLI contracts; no catalogue access."""
import csv
from dataclasses import replace
from pathlib import Path
import socket
import subprocess
import sys

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import numpy as np
import pytest

import barnard_motion as b

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def source():
    return b.load_source()


def test_catalogue_epoch_and_formal_covariance(source):
    result = b.predict(source, source.epoch)
    initial = SkyCoord(source.values[0] * u.deg, source.values[1] * u.deg)
    assert initial.separation(result.coordinate).to_value(u.mas) < 0.001
    np.testing.assert_allclose(np.sqrt(np.diag(result.covariance_mas2)),
                               np.sqrt(np.diag(source.covariance))[:2], rtol=0.01)
    np.testing.assert_allclose(result.covariance_mas2, source.covariance[:2, :2], rtol=0.01)
    np.testing.assert_allclose(result.jacobian[:, :2], np.eye(2), atol=1e-6)


@pytest.mark.parametrize("declination", [4.739420051112412, 60.0])
def test_short_motion_and_cosine(source, declination):
    values = source.values.copy()
    values[1] = declination
    local = replace(source, values=values)
    dt = 0.01
    initial = b.coordinate_at(local, local.epoch)
    moved = b.coordinate_at(local, local.epoch + dt * u.yr)
    np.testing.assert_allclose(b.tangent_offsets_mas(initial, moved),
                               values[3:5] * dt, atol=0.001, rtol=0)
    raw_ra_mas = (moved.ra - initial.ra).to_value(u.mas)
    assert raw_ra_mas == pytest.approx(values[3] * dt / np.cos(np.deg2rad(declination)), abs=0.001)


def test_short_covariance_has_position_pm_correlations(source):
    dt = 0.01
    transform = np.array([[1, 0, 0, dt, 0, 0], [0, 1, 0, 0, dt, 0]])
    expected = transform @ source.covariance @ transform.T
    actual = b.predict(source, source.epoch + dt * u.yr).covariance_mas2
    np.testing.assert_allclose(actual, expected, rtol=1e-4, atol=1e-10)


def test_rv_and_perspective_acceleration_are_used(source):
    time = b.parse_time("2025-10-01")
    rv_zero = source.values.copy()
    rv_zero[5] = 0
    difference = b.coordinate_at(source, time).separation(
        b.coordinate_at(replace(source, values=rv_zero), time)).to_value(u.mas)
    assert difference > 10  # Barnard's approaching velocity matters by many mas.
    prediction = b.predict(source, time)
    assert np.linalg.norm(prediction.jacobian[:, 5]) > 0.01
    no_rv_error = source.covariance.copy()
    no_rv_error[5, 5] = 0
    assert np.trace(prediction.covariance_mas2) > np.trace(
        b.predict(replace(source, covariance=no_rv_error), time).covariance_mas2)


@pytest.mark.parametrize("date", ["2025-04-01", "2025-10-01"])
def test_numerical_jacobian_step_stability(source, date):
    time = b.parse_time(date)
    standard = b.predict(source, time).covariance_mas2
    for multiplier in [0.5, 2.0]:
        np.testing.assert_allclose(b.predict(source, time, step_scale=multiplier).covariance_mas2,
                                   standard, rtol=1e-4, atol=1e-9)


def test_same_instant_in_different_formats(source):
    time = b.parse_time("2025-04-01T00:00:00Z")
    representations = [b.parse_time(f"{time.mjd:.14f}", "mjd", "utc"),
                       b.parse_time(f"{time.tcb.jyear:.14f}", "jyear", "tcb"),
                       b.parse_time(time.tcb.isot, "isot", "tcb")]
    expected = b.predict(source, time)
    for other in representations:
        assert abs((other.tai - time.tai).to_value(u.s)) < 1e-4
        actual = b.predict(source, other)
        assert expected.coordinate.separation(actual.coordinate).to_value(u.mas) < 0.001
        np.testing.assert_allclose(actual.covariance_mas2, expected.covariance_mas2, rtol=1e-4)


def test_midpoint_uses_continuous_time_through_leap_second():
    row = {"obs_id": "00012345678901234567890", "time_format": "iso", "time_scale": "utc",
           "t_min": "2016-12-31T23:59:59", "t_max": "2017-01-01T00:00:00"}
    obs = b.observation_from_row(row)
    assert obs.time.utc.isot == "2016-12-31T23:59:60.000000000"
    assert obs.obs_id == row["obs_id"]
    mjd = b.observation_from_row({"obs_id": "interval", "time_format": "mjd", "time_scale": "tcb",
                                 "t_min": "60800.1", "t_max": "60800.3"})
    assert mjd.time.tcb.mjd == pytest.approx(60800.2, abs=1e-10)


@pytest.mark.parametrize("value,fmt,scale,explanation", [
    ("2025-02-30", None, None, "Некорректная дата"),
    ("2025-13-01", None, None, "Некорректная дата"),
    ("2025-01-01T24:00:00", None, None, "ISO-дата"),
    ("2025-01-01T12:00:00+05:00", None, None, "ISO-дата"),
    ("60800", None, None, "числовые даты"),
    (60800, "mjd", None, "шкалу времени явно"),
    ("nan", "mjd", "utc", "конечное число"),
    ("inf", "jyear", "tcb", "конечное число"),
    ("2025-04-01", "iso", "local", "шкала"),
    ("2025-04-01", "iso", "ut1", "UT1"),
    ("2025-04-01Z", "iso", "tcb", "Z"),
    ("2025", "year", "utc", "формат"),
])
def test_invalid_times_are_explained(value, fmt, scale, explanation):
    with pytest.raises(ValueError, match=explanation):
        b.parse_time(value, fmt, scale)


@pytest.mark.parametrize("updates,explanation", [
    ({"time": "2025-04-01", "t_min": "2025-04-01", "t_max": "2025-04-02"}, "Неоднозначная"),
    ({"t_min": "2025-04-02", "t_max": "2025-04-01"}, "Обратный интервал"),
    ({"t_min": "2025-04-01"}, "полная пара"),
    ({}, "Требуется time"),
    ({"time": "2025-04-01", "time_scale": ""}, "time_scale"),
    ({"time": "2025-04-01", "obs_id": ""}, "obs_id"),
    ({"time": "2025-04-01", "is_demo": "maybe"}, "is_demo"),
])
def test_invalid_observations_are_explained(updates, explanation):
    with pytest.raises(ValueError, match=explanation):
        b.observation_from_row({"obs_id": "test-only", "time_format": "iso", "time_scale": "utc", **updates})


@pytest.mark.parametrize("key,value", [("radial_velocity_km_s", ""), ("ra_dec_corr", ""),
                                        ("parallax_mas", "-1"), ("ra_error_basis", "sigma_alpha"),
                                        ("ref_epoch_scale", "UTC"), ("pmra_error_mas_yr", "0"),
                                        ("source_id", "4472832130942576000")])
def test_missing_or_invalid_source_parameters(source, tmp_path, key, value):
    row = {**source.row, key: value}
    path = tmp_path / "catalog.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError):
        b.load_source(path)


def test_bad_csv_structure_and_row_context(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("obs_id,time_format,time_scale,time\nx,iso,utc,2025-02-30\n")
    with pytest.raises(ValueError, match="строка 2.*x.*Некорректная дата"):
        b.load_observations(path)
    path.write_text("obs_id,time_format,time_scale,time\nx,iso,utc\n")
    with pytest.raises(ValueError, match="число полей"):
        b.load_observations(path)
    path.write_text("obs_id,time_format\n")
    with pytest.raises(ValueError, match="time_scale"):
        b.load_observations(path)


def test_single_batch_cli_and_exact_identifiers(source, tmp_path):
    path = tmp_path / "observations.csv"
    identifier = "000123456789012345678901234567890"
    path.write_text("obs_id,time_format,time_scale,time,detector,data_release,access_url\n"
                    f"{identifier},iso,utc,2025-04-01T00:00:00,D6,test-release,https://example.invalid/test\n")
    single, batch = tmp_path / "single.csv", tmp_path / "batch.csv"
    for args in [["--date", "2025-04-01T00:00:00", "--output", str(single)],
                 ["--observations", str(path), "--output", str(batch)]]:
        run = subprocess.run([sys.executable, str(ROOT / "barnard_motion.py"), *args],
                             cwd=tmp_path, text=True, capture_output=True)
        assert run.returncode == 0, run.stderr
    with single.open() as f, batch.open() as g:
        one, many = next(csv.DictReader(f)), next(csv.DictReader(g))
    assert one["obs_id"] == ""
    assert many["obs_id"] == identifier
    assert many["source_id"] == b.SOURCE_ID
    assert many["detector"] == "D6" and many["data_release"] == "test-release"
    assert many["parallax_applied"] == "false"
    for key in b.OUTPUT_COLUMNS:
        if key not in {"obs_id", "detector", "data_release", "access_url"}:
            assert one[key] == many[key], key
    expected = b.calculate_observation(source, b.Observation(b.parse_time("2025-04-01")))
    assert float(one["ra_deg"]) == expected["ra_deg"]


def test_offline_calculation(source, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Calculation attempted network access")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "getaddrinfo", forbidden)
    rows = b.calculate_observations(source, b.load_observations(ROOT / "data/barnard_observations_demo.csv"))
    assert len(rows) == 2
    assert all(row["is_demo"] == "true" and row["obs_id"] == "" for row in rows)
    assert all(row["position_type"] == "barycentric_astrometric" for row in rows)
