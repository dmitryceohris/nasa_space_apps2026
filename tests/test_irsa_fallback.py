"""TAP outage recovery through official SIA2, without stale results or network access."""
import io
import json
from types import SimpleNamespace

import pytest
from astropy.io.votable import from_table, writeto
from astropy.table import Table

from api import irsa
from api.errors import ApiError


def sia_document(status="OK", empty=False):
    table = Table(names=["access_url", "obs_id", "obs_collection", "energy_bandpassname",
                         "t_min", "t_max", "em_min", "em_max"],
                  dtype=["U160", "U24", "U24", "U16", float, float, float, float])
    if not empty:
        for collection, detector in [("spherex_qr3", 3), ("spherex_dr1", 3), ("spherex_dr1", 2)]:
            table.add_row((f"https://irsa.ipac.caltech.edu/ibe/data/{collection}/level2_test.fits",
                           "same-observation", collection, f"SPHEREx-D{detector}", 61000., 61000.2, 1.6e-6, 2.4e-6))
    votable = from_table(table)
    from astropy.io.votable.tree import Info
    votable.resources[0].infos.append(Info(name="QUERY_STATUS", value=status))
    stream = io.BytesIO()
    writeto(votable, stream)
    return stream.getvalue()


@pytest.fixture(autouse=True)
def isolated_archive(monkeypatch, tmp_path):
    monkeypatch.setattr(irsa, "CACHE", tmp_path)
    monkeypatch.setattr(irsa, "_tap_retry_at", 0.)


def test_tap_failure_uses_sia_and_cache_preserves_provenance(monkeypatch):
    calls = []

    def send(method, url, **kwargs):
        calls.append((method, url, kwargs))
        if url.endswith("/TAP/sync"):
            assert method == "GET" and kwargs["retries"] == 1 and kwargs["timeout"] == 12
            raise ApiError(502, "upstream", "TAP unavailable")
        assert url.endswith("/SIA") and method == "GET"
        assert kwargs["params"]["COLLECTION"] == ["spherex_dr1", "spherex_qr3"]
        return SimpleNamespace(content=sia_document())

    monkeypatch.setattr(irsa, "_send", send)
    rows = irsa._tap_exposures(10., 20., 3)
    assert rows.source == "sia2" and rows.warnings == [irsa.SIA_LIMITATION]
    assert len(rows) == 1 and "spherex_dr1" in rows[0]["uri"]
    assert rows[0]["mjd"] == pytest.approx(61000.1)
    assert rows[0]["lam"] == pytest.approx([1.6, 2.4])
    cached = irsa._tap_exposures(10., 20., 3)
    assert cached == rows and cached.source == "sia2" and cached.warnings == rows.warnings
    assert len(calls) == 2
    irsa._tap_exposures(11., 20., 3)
    assert len(calls) == 3 and calls[-1][1].endswith("/SIA")  # TAP cooling down


def test_tap_is_probed_again_after_cooldown(monkeypatch):
    monkeypatch.setattr(irsa, "_tap_retry_at", 100.)
    monkeypatch.setattr(irsa.time, "monotonic", lambda: 101.)
    calls = []

    def send(method, url, **kwargs):
        calls.append(url)
        return SimpleNamespace(text="uri,mjd,lo,hi\nibe/level2_test.fits,61000,1e-6,2e-6\n")

    monkeypatch.setattr(irsa, "_send", send)
    rows = irsa._tap_exposures(10., 20., 3)
    assert rows.source == "tap" and rows.warnings == [] and len(rows) == 1
    assert len(calls) == 1 and calls[0].endswith("/TAP/sync")


@pytest.mark.parametrize("status", ["ERROR", "OVERFLOW"])
def test_sia_error_or_truncation_is_not_mistaken_for_no_coverage(monkeypatch, status):
    monkeypatch.setattr(irsa, "_tap_retry_at", float("inf"))
    monkeypatch.setattr(irsa, "_send", lambda *a, **kw: SimpleNamespace(content=sia_document(status)))
    with pytest.raises(ApiError) as error:
        irsa._tap_exposures(10., 20., 3)
    assert error.value.status == 502
    assert list(irsa.CACHE.glob("exp_*.json")) == []


def test_empty_sia_result_keeps_its_source_and_warning(monkeypatch):
    monkeypatch.setattr(irsa, "_tap_retry_at", float("inf"))
    monkeypatch.setattr(irsa, "_send", lambda *a, **kw: SimpleNamespace(content=sia_document(empty=True)))
    rows = irsa._tap_exposures(10., 20., 3)
    assert rows == [] and rows.source == "sia2" and rows.warnings
    saved = json.loads(next(irsa.CACHE.glob("exp_*.json")).read_text())
    assert saved["source"] == "sia2" and saved["rows"] == []


def test_direct_sia_fits_url_does_not_get_parsed_as_datalink(monkeypatch):
    monkeypatch.setattr(irsa, "_send", lambda *a, **kw: pytest.fail("direct FITS is not a datalink document"))
    url = "https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/level2_test.fits"
    assert irsa._fits_url(url) == url
