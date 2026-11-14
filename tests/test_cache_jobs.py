"""Bounded-job and cache freshness checks with controlled clocks and no network."""
import io
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits

from api import cache, irsa, jobs
from api.errors import ApiError


@pytest.mark.parametrize('version', [8, 9, 10])
def test_previous_result_schemas_are_not_served(monkeypatch, tmp_path, version):
    monkeypatch.setattr(cache, 'CACHE', tmp_path)
    body = {'status': 'ready', 'frames': [{}]}
    cache.write('old', body)
    path = tmp_path / 'old' / 'result.json'
    result = json.loads(path.read_text())
    result['v'] = version
    path.write_text(json.dumps(result))
    assert cache.read('old') is None
    cache.write('current', body)
    assert cache.read('current') == body


def test_expired_partial_and_empty_results_refresh_sooner(monkeypatch, tmp_path):
    monkeypatch.setattr(cache, 'CACHE', tmp_path)
    monkeypatch.setattr(cache.time, 'time', lambda: 1000)
    for key, body in [('complete', {'status':'ready','frames':[{}]}),
                      ('partial', {'status':'ready','data_status':'partial','frames':[{}]}),
                      ('empty', {'status':'ready','frames':[]})]:
        cache.write(key, body)
        os.utime(tmp_path/key/'result.json', (900, 900))
    assert cache.read('complete') is not None
    assert cache.read('partial') is None and cache.read('empty') is None
    os.utime(tmp_path/'complete'/'result.json', (-3000, -3000))
    assert cache.read('complete') is None


def test_completed_worker_during_cache_miss_does_not_start_duplicate(monkeypatch):
    reads = iter([None, {'status':'ready','frames':[{}]}])
    monkeypatch.setattr(cache, 'read', lambda key: next(reads))
    monkeypatch.setattr(jobs, 'JOBS', {})
    monkeypatch.setattr(jobs.POOL, 'submit', lambda *args: pytest.fail('duplicate worker'))
    assert jobs._poll('key', lambda job: None)[0] == 200


def test_same_query_shares_one_queued_job_and_queue_is_bounded(monkeypatch):
    monkeypatch.setattr(cache, 'read', lambda key: None)
    monkeypatch.setattr(jobs, 'JOBS', {})
    monkeypatch.setenv('MAX_PENDING_JOBS', '1')
    submitted = []
    monkeypatch.setattr(jobs.POOL, 'submit', lambda *args: submitted.append(args))
    first = jobs._poll('one', lambda job: None)
    second = jobs._poll('one', lambda job: None)
    assert first[1]['job_id'] == second[1]['job_id'] and len(submitted) == 1
    assert first[1]['stage'] == 'queued'
    with pytest.raises(ApiError) as error:
        jobs._poll('two', lambda job: None)
    assert error.value.code == 'busy' and error.value.status == 503


def test_expired_queued_job_never_runs_or_publishes(monkeypatch):
    monkeypatch.setattr(jobs.time, 'monotonic', lambda: 100)
    monkeypatch.setenv('MAX_JOB_SECONDS', '10')
    job = jobs.Job('key', 'id', created=0)
    jobs._run(job, lambda job: pytest.fail('expired work started'))
    assert job.status == 'error' and 'время' in job.message


def test_exposure_cache_expires_and_tables_are_explicit(monkeypatch, tmp_path):
    monkeypatch.setattr(irsa, 'CACHE', tmp_path)
    monkeypatch.setattr(irsa.time, 'time', lambda: 1000)
    path = tmp_path/'exp_10.000_20.000_d3.json'
    path.write_text('[]')
    os.utime(path, (900, 900))
    requested = []
    def send(*args, **kwargs):
        requested.append(args)
        return SimpleNamespace(text='uri,mjd,lo,hi\nibe/level2_example.fits,61000,1e-6,2e-6\n')
    monkeypatch.setattr(irsa, '_send', send)
    assert len(irsa.exposures(10,20,3)) == 1 and len(requested) == 1
    monkeypatch.setenv('SPHEREX_SOURCE', 'tables')
    monkeypatch.setattr(irsa, '_table_rows', lambda: [])
    assert irsa.exposures(10,20,3) == [] and len(requested) == 1


def test_broken_fits_cache_is_replaced_atomically(monkeypatch, tmp_path):
    monkeypatch.setattr(irsa, 'CACHE', tmp_path)
    # No network: a valid real FITS container replaces an interrupted file.
    hdul = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.ones((5,5)), name='IMAGE')])
    stream = io.BytesIO(); hdul.writeto(stream)
    monkeypatch.setattr(irsa, '_send', lambda *a, **kw: SimpleNamespace(content=stream.getvalue()))
    result = irsa.cutout('ibe/level2_example.fits', 10,20,.1)
    result.write_bytes(b'incomplete')
    assert irsa.cutout('ibe/level2_example.fits', 10,20,.1) == result
    assert fits.getdata(result, 'IMAGE').shape == (5,5)
    assert list(tmp_path.glob('tmp*')) == []


def test_timed_out_job_remains_reserved_while_inflight_calls_finish(monkeypatch):
    job = jobs.Job('key', 'id', created=0, started=True)
    monkeypatch.setattr(jobs, 'JOBS', {'key':job})
    monkeypatch.setattr(cache, 'read', lambda key: None)
    monkeypatch.setattr(jobs.time, 'monotonic', lambda:100)
    monkeypatch.setenv('MAX_JOB_SECONDS','10')
    monkeypatch.setattr(jobs.POOL, 'submit', lambda *args:pytest.fail('duplicated in-flight work'))
    with pytest.raises(ApiError) as error:
        jobs._poll('key', lambda job:None)
    assert error.value.code == 'timeout' and jobs.JOBS['key'] is job
