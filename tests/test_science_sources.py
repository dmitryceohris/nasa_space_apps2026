"""Check catalogue identity/epoch against saved primary-source responses, without network."""
import csv
import hashlib
import json
from pathlib import Path

import pytest

from api import resolve

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "data/validation/science-acceptance"


def test_primary_source_snapshots_and_saved_observations_are_intact():
    for source in json.loads((SOURCES / "sources.json").read_text()):
        assert hashlib.sha256((SOURCES / source["file"]).read_bytes()).hexdigest() == source["sha256"]
    evidence = json.loads((SOURCES / "evidence.json").read_text())
    assert len(evidence["frames"]) == 43
    assert len(evidence["wise_exact_coordinate_frames"]) == 4
    for row in evidence["frames"] + evidence["wise_exact_coordinate_frames"]:
        assert hashlib.sha256((ROOT / row["response"]).read_bytes()).hexdigest() == row["response_sha256"]


def test_five_catalogue_rows_match_primary_parameters_and_keep_component_caveats():
    gaia = {r["source_id"]: r for r in csv.DictReader((SOURCES / "gaia.csv").open())}
    parameters = list(csv.DictReader((ROOT / "OBJECT_PARAMETERS.csv").open()))
    mapping = {"ra_deg": "ra", "dec_deg": "dec", "pmra_mas_yr": "pmra", "pmdec_mas_yr": "pmdec",
               "parallax_mas": "parallax", "ra_error_mas": "ra_error", "dec_error_mas": "dec_error",
               "pmra_error_mas_yr": "pmra_error", "pmdec_error_mas_yr": "pmdec_error",
               "parallax_error_mas": "parallax_error", "radial_velocity_km_s": "radial_velocity"}
    for row in parameters:
        if row["source_id"] in gaia:
            primary = gaia[row["source_id"]]
            for ours, theirs in mapping.items():
                if row[ours]:
                    assert float(row[ours]) == pytest.approx(float(primary[theirs]), rel=1e-12)
                else:
                    assert not primary[theirs]
            assert row["ref_epoch"] == primary["ref_epoch"] and row["ref_epoch_format"] == "jyear"
    lines = (SOURCES / "wise-plx.tsv").read_text().splitlines()
    header = next(line for line in lines if line.startswith("recno\t")).split("\t")
    primary = dict(zip(header, next(line for line in lines if "\t0855-0714\t" in line).split("\t")))
    wise = next(r for r in parameters if r["source_id"] == "0855-0714")
    for ours, theirs in {"ra_deg": "RA_ICRS", "dec_deg": "DE_ICRS", "ref_epoch": "MJD",
                        "pmra_mas_yr": "pmRA", "pmdec_mas_yr": "pmDE", "parallax_mas": "plx"}.items():
        assert float(wise[ours]) == float(primary[theirs])
    assert wise["ref_epoch_format"] == "mjd" and wise["ref_epoch_scale"] == ""
    eps = resolve.resolve("epsindib")
    assert eps["source_id"] == "6412596012146801152" and "not the Ba/Bb barycentre" in next(
        r["astrometric_target"] for r in parameters if r["source_id"] == eps["source_id"])
    rw = resolve.resolve("rwaur")
    assert rw["source_id"] == "156430817820015232" and "eastern" in rw["component"]
    assert "Gaia DR3 156430817820015232" in (SOURCES / "rw-B-simbad.txt").read_text()
    assert "Gaia DR3 156430822114424576" in (SOURCES / "rw-A-simbad.txt").read_text()


@pytest.mark.parametrize("name,filename,canonical,alias", [
    ("Proxima Cen", "proxima-cen", "NAME Proxima Centauri", "Gaia DR3 5853498713190525696"),
    ("61 Cyg A", "61-cyg-a", "*  61 Cyg A", "Gaia DR3 1872046609345556480"),
    ("Vega", "vega", "* alf Lyr", "HD 172167"),
])
def test_saved_sesame_identity_does_not_invent_solution_epoch(name, filename, canonical, alias):
    text = (ROOT / f"data/validation/responses/{filename}-sesame.txt").read_text()
    result = resolve._parse(name, [text])
    assert result["source_id"] == canonical and alias in text
    assert result["ref_epoch"] is result["ref_epoch_format"] is result["ref_epoch_scale"] is None
    assert "reference_epoch_unverified" in result["assumptions"]
    assert result["pmra_masyr"] is not None and result["pmdec_masyr"] is not None
    assert resolve.propagate(result, 61358) == (result["ra"], result["dec"])
