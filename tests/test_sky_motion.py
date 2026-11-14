"""Physical invariants and small local quality/selection checks; no network."""
import csv
from dataclasses import replace
import socket

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, SkyCoord
from astropy.io import fits
from astropy.time import Time
import numpy as np
import pytest

from barnard_motion import load_observations
import sky_motion as sky


@pytest.fixture
def targets():
    return sky.load_targets()


@pytest.mark.parametrize("name", sky.TARGET_IDS.values())
def test_five_reference_positions_and_assumptions(targets, name):
    target = targets[name]
    prediction = sky.predict_position(target, target.epoch)
    reference = SkyCoord(target.values[0] * u.deg, target.values[1] * u.deg)
    assert sky.angular_position(prediction.barycentric).separation(reference).mas < 0.001
    assert np.isfinite(prediction.geocentric.ra.deg)
    assert np.isfinite(prediction.geocentric.dec.deg)
    assert ("rv_assumed_zero" in target.assumptions) == (name != "Barnard's Star")
    if name == "WISE 0855-0714":
        assert "reference_epoch_scale_assumed_utc" in target.assumptions
        assert target.epoch.utc.mjd == pytest.approx(57633.39)
        assert target.row["source_id"] == "0855-0714"


def test_pmra_has_cos_declination(targets):
    target = replace(targets["Barnard's Star"], values=np.array([30., 60., 10., 1000., 0., 0.]))
    position = sky.predict_position(target, target.epoch + 0.1 * u.yr).barycentric
    assert (position.ra.deg - 30) * 3600 == pytest.approx(0.2, abs=1e-5)


@pytest.mark.parametrize("axis,component", [(1, "east"), (2, "north")])
def test_parallax_sign_and_one_au_scale(targets, monkeypatch, axis, component):
    target = replace(targets["Barnard's Star"], values=np.array([0., 0., 1000., 0., 0., 0.]))
    vector = np.zeros(3)
    vector[axis] = 1
    monkeypatch.setattr(sky, "get_body_barycentric", lambda *args, **kwargs: CartesianRepresentation(vector * u.au))
    prediction = sky.predict_position(target, target.epoch)
    east, north = sky.angular_position(prediction.barycentric).spherical_offsets_to(sky.angular_position(prediction.geocentric))
    assert {"east": east.arcsec, "north": north.arcsec}[component] == pytest.approx(-1.0, abs=1e-7)


def test_original_barnard_regression(targets):
    observations = load_observations(sky.ROOT / "data/irsa_barnard_20261114/observations.csv")
    positions = [sky.predict_position(targets["Barnard's Star"], obs.time) for obs in observations]
    bary = [sky.angular_position(position.barycentric) for position in positions]
    geo = [sky.angular_position(position.geocentric) for position in positions]
    assert bary[0].separation(bary[-1]).arcsec == pytest.approx(5.594, abs=0.002)
    assert geo[0].separation(geo[-1]).arcsec == pytest.approx(5.455, abs=0.002)
    assert bary[0].spherical_offsets_to(bary[-1])[0].arcsec < 0
    assert geo[0].spherical_offsets_to(geo[-1])[0].arcsec > 0


def test_epoch_clusters_and_continuous_windows():
    groups = sky.select_epochs([{"mjd": value} for value in [100, 101, 135, 136, 200, 235]])
    assert len(groups) == 3
    assert [[row["mjd"] for row in group] for _, group in groups] == [[100, 101], [200], [235]]
    continuous = sky.select_epochs([{"mjd": value} for value in [100, 105, 110, 115, 120, 125]])
    assert [name for name, _ in continuous] == ["window_1", "window_2", "window_3"]
    assert [len(group) for _, group in continuous] == [2, 2, 2]
    assert len(sky.select_epochs([{"mjd": 100}, {"mjd": 130}, {"mjd": 161}])) == 2


def test_spherical_footprint_and_ra_wrap():
    region = "POLYGON ICRS 359 -1 1 -1 1 1 359 1 359 -1"
    assert sky.footprint_margin(region, SkyCoord(0 * u.deg, 0 * u.deg)) > 3599
    assert sky.footprint_margin(region, SkyCoord(2 * u.deg, 0 * u.deg)) < 0


@pytest.mark.parametrize("bad_name", ["OVERFLOW", "NONLINEAR", "NONFUNC", "MISSING_DATA"])
def test_local_flags_use_product_bit_definitions(bad_name):
    # Deliberately use arbitrary bits, rather than the current archive bit numbers.
    definitions = {"MSKN0003": "OVERFLOW", "MSKN0007": "NONLINEAR", "MSKN0009": "NONFUNC", "MSKN0011": "MISSING_DATA",
                   "MSKN0000": "SUR_ERROR", "MSKN0002": "HOT", "MSKN0005": "COLD"}
    header = fits.Header(definitions)
    data = np.zeros((15, 15), dtype=np.uint32)
    data[0, 0] = 1 << 3  # A bad pixel outside the target region must not reject it.
    assert sky.assess_flags(data, header, 7, 7)[0] == "usable"
    bit = int(next(key[4:] for key, value in definitions.items() if value == bad_name))
    data[6, 8] = 1 << bit
    status, evidence = sky.assess_flags(data, header, 7, 7)
    assert status == "flagged"
    assert evidence["bad_flag_counts"][bad_name] == 1
    assert sky.assess_flags(data, header, 0, 0)[0] == "unchecked"
    assert sky.assess_flags(data, fits.Header(), 7, 7)[0] == "unchecked"


def test_time_precedence_scale_and_fallback():
    row = {"t_min": "60914", "t_max": "60916"}
    header = fits.Header({"MJD-AVG": 60915.25, "DATE-AVG": "2025-08-27", "TIMESYS": "UTC"})
    moment, method, assumptions = sky.observation_time(header, row)
    assert moment.utc.mjd == pytest.approx(60915.25)
    assert method == "fits_mjd_avg" and not assumptions
    moment, method, assumptions = sky.observation_time(fits.Header(), row)
    assert moment.utc.mjd == pytest.approx(60915)
    assert method == "sia_midpoint"
    assert "sia_time_scale_assumed_utc" in assumptions
    with pytest.raises(ValueError, match="Reversed"):
        sky.observation_time(fits.Header(), {"t_min": "60916", "t_max": "60914"})


def test_offline_batch_and_csv(targets, monkeypatch, tmp_path):
    def forbid_network(*args, **kwargs):
        raise AssertionError("Offline computation attempted a network connection")
    monkeypatch.setattr(socket.socket, "connect", forbid_network)
    rows = [{"target": name, "mjd": "60915", "quality_status": "unchecked"} for name in targets]
    rows.append({"target": "Barnard's Star", "mjd": "61110", "quality_status": "flagged"})
    predictions = sky.calculate_predictions(rows, targets)
    assert len(predictions) == 6
    first_rows = [row for row in predictions if row["mjd"] == "60915"]
    assert all(abs(row["shift_px"]) < 1e-8 for row in first_rows)
    last = next(row for row in predictions if row["mjd"] == "61110")
    assert last["shift_px"] > 0.8
    assert last["quality_status"] == "flagged"
    output = tmp_path / "predictions.csv"
    sky.write_csv(predictions, output, sky.PREDICTION_COLUMNS)
    with output.open(newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == list(sky.PREDICTION_COLUMNS)
        assert len(list(reader)) == 6


def test_unavailable_observation_is_not_invented(targets):
    result = sky.calculate_predictions([{"target": "R Dra", "quality_status": "unchecked", "assumptions": "no_observations_prepared"}], targets)[0]
    assert result["ra_geo"] == "" and result["mjd"] == ""
    assert result["quality_status"] == "unchecked"


def test_known_bad_candidate_is_not_hidden_by_unchecked_replacement():
    attempts = [{"result": {"quality_status": status}} for status in ["unchecked", "flagged"]]
    assert sky.choose_attempt(attempts)["quality_status"] == "flagged"
    attempts.append({"result": {"quality_status": "usable"}})
    assert sky.choose_attempt(attempts)["quality_status"] == "usable"
