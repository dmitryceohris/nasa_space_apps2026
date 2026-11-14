"""The spectrum/library archive path must retry upstream failures and bound waits."""
from types import SimpleNamespace

import pytest
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table
from pyvo.dal.exceptions import DALServiceError

from Exposure import irsa
from api import jobs
from api.errors import ApiError


@pytest.mark.parametrize("failure", [DALServiceError("proxy", code=502),
                                    DALServiceError("busy", code=503), requests.Timeout("read")])
def test_temporary_sia_failure_retries_before_returning_data(monkeypatch, failure):
    table = Table({"obs_id": ["real-shape"]})
    calls, pauses = [], []

    def search(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            raise failure
        return SimpleNamespace(to_table=lambda: table)

    monkeypatch.setattr(irsa, "SIA2Service", lambda *a, **kw: SimpleNamespace(search=search))
    monkeypatch.setattr(irsa.time, "sleep", pauses.append)
    assert irsa._query_sia(SkyCoord(10, 20, unit="deg"), irsa.DEFAULT_COLLECTIONS) is table
    assert len(calls) == 2 and pauses == [2]


@pytest.mark.parametrize("status, attempts", [(503, 3), (403, 1)])
def test_failed_sia_search_is_reported_as_archive_failure(monkeypatch, status, attempts):
    calls = []

    def search(**kwargs):
        calls.append(kwargs)
        raise DALServiceError("upstream failure", code=status)

    monkeypatch.setattr(irsa, "SIA2Service", lambda *a, **kw: SimpleNamespace(search=search))
    monkeypatch.setattr(irsa.time, "sleep", lambda seconds: None)
    with pytest.raises(irsa.ArchiveQueryError):
        irsa._query_sia(SkyCoord(10, 20, unit="deg"), irsa.DEFAULT_COLLECTIONS)
    assert len(calls) == attempts


def test_pyvo_session_does_not_leave_waits_unbounded(monkeypatch):
    waits = []

    def request(self, method, url, **kwargs):
        waits.append(kwargs.get("timeout"))
        return object()

    monkeypatch.setattr(requests.Session, "request", request)
    with irsa._ArchiveSession() as session:
        session.get("https://irsa.example/SIA")
        session.get("https://irsa.example/SIA", timeout=None)
        session.get("https://irsa.example/SIA", timeout=.25)
    assert all(isinstance(timeout, (float, int)) and 0 < timeout <= 45 for timeout in waits)
    assert waits[-1] == .25


def test_spectrum_maps_archive_failure_to_upstream_error(monkeypatch):
    def fail(*args, **kwargs):
        raise irsa.ArchiveQueryError("IRSA SIA2 unavailable")

    monkeypatch.setattr(jobs, "build_spectrum", fail)
    with pytest.raises(ApiError) as error:
        jobs._spectrum_job(jobs.Job("test", "id"), {"center": [10., 20.]})
    assert error.value.status == 502 and error.value.code == "upstream"
