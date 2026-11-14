"""Public errors and bounded throttling; no external requests."""
import os
os.environ['PREWARM'] = '0'
import pytest
from fastapi.testclient import TestClient
from api import main, jobs
from api.errors import ApiError, PUBLIC_INTERNAL_MESSAGE
from api.limits import RequestLimiter


@pytest.fixture(autouse=True)
def isolated_limiter(monkeypatch):
    monkeypatch.setattr(main, 'request_limiter', RequestLimiter(10000))


def test_unexpected_error_is_private_but_logged(monkeypatch, caplog):
    def crash(*args):
        raise RuntimeError('PRIVATE-secret-path')
    monkeypatch.setattr(jobs, 'view', crash)
    response = TestClient(main.app, raise_server_exceptions=False).get('/api/view?ra=10&dec=20')
    assert response.status_code == 500
    assert response.json()['error']['message'] == PUBLIC_INTERNAL_MESSAGE
    assert 'PRIVATE-secret-path' not in response.text
    assert 'PRIVATE-secret-path' in caplog.text


def test_background_exception_is_private_but_logged(monkeypatch, caplog):
    def crash(job):
        raise RuntimeError('PRIVATE-background-path')
    job = jobs.Job('test', 'test')
    jobs._run(job, crash)
    assert job.status == 'error' and job.message == PUBLIC_INTERNAL_MESSAGE
    assert 'PRIVATE-background-path' in caplog.text


def test_throttled_requests_do_not_start_work_and_refill(monkeypatch):
    now = [0.0]
    monkeypatch.setattr(main, 'request_limiter', RequestLimiter(2, clock=lambda:now[0]))
    calls = []
    monkeypatch.setattr(jobs, 'view', lambda p:(calls.append(p) or (200, {'status':'ready','frames':[]})))
    client = TestClient(main.app)
    for _ in range(2):
        assert client.get('/api/view?ra=10&dec=20').status_code == 200
    response = client.get('/api/view?ra=10&dec=20', headers={'X-Forwarded-For':'fake-new-peer','Origin':'https://example.org'})
    assert response.status_code == 429 and response.headers['retry-after'] == '30'
    assert response.json()['error']['code'] == 'rate_limited' and len(calls) == 2
    assert response.headers['access-control-allow-origin'] == '*'
    assert 'Retry-After' in response.headers['access-control-expose-headers']
    assert client.get('/api/health').status_code == 200
    now[0] = 30
    assert client.get('/api/view?ra=10&dec=20').status_code == 200


def test_buckets_are_per_peer_and_memory_is_bounded():
    limiter = RequestLimiter(1, max_clients=2)
    assert limiter.retry_after('a') is None
    assert limiter.retry_after('a') is not None
    assert limiter.retry_after('b') is None
    limiter.retry_after('c')
    assert len(limiter.clients) == 2


def test_busy_response_preserves_error_shape_and_retry_header(monkeypatch):
    def busy(*args):
        raise ApiError(503, 'busy', 'Queue full')
    monkeypatch.setattr(jobs, 'view', busy)
    response = TestClient(main.app).get('/api/view?ra=10&dec=20')
    assert response.status_code == 503 and int(response.headers['retry-after']) > 0
    assert response.json()['error'] == {'code':'busy','message':'Queue full'}
