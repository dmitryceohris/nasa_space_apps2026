"""Offline resolver contract and common-model integration checks."""
from urllib.parse import unquote, urlsplit

import httpx
import pytest
from astropy.time import Time

import sky_motion
from api import resolve


@pytest.mark.parametrize('name', ["Cl* NGC 123+45 #2", "Barnard's other star"])
def test_complete_name_is_encoded(monkeypatch, name):
    urls = []
    def get(url, **kwargs):
        urls.append(url)
        return httpx.Response(200, text='Nothing found')
    monkeypatch.setattr(resolve.httpx, 'get', get)
    resolve._fetch_sesame(name)
    parts = urlsplit(urls[0])
    assert not parts.fragment
    assert unquote(parts.query) == name
    assert '+' not in parts.query


@pytest.mark.parametrize('tid', ['barnard', 'epsindib', 'wise0855', 'rwaur', 'rdra'])
def test_catalog_provenance_unknowns_and_common_prediction(monkeypatch, tid):
    monkeypatch.setattr(resolve, '_fetch_sesame', lambda _: pytest.fail('catalog queried Sesame'))
    monkeypatch.setitem(resolve._mem, tid, {'ra': 0, 'dec': 0})
    record = resolve.resolve(tid)
    assert record['id'] == tid
    assert isinstance(record['source_id'], str)
    assert record['source_url'] and record['ref_epoch_format']
    target = sky_motion.load_targets()[record['name']]
    assert record['ref_epoch'] == float(target.row['ref_epoch'])
    assert record['rv_kms'] == (float(target.row['radial_velocity_km_s']) if target.row['radial_velocity_km_s'] else None)
    for mjd in (60000, 60903, 61200):
        predicted = sky_motion.predict_position(target, Time(mjd, format='mjd', scale='utc')).geocentric
        assert resolve.propagate(record, mjd) == pytest.approx((predicted.ra.deg, predicted.dec.deg), abs=1e-10)


def test_sesame_does_not_invent_epoch_or_mix_catalogs():
    text = '#=Sc=Simbad\n%J 10 -20\n%P 1 2\n#=V=Vizier\n%J 30 40\n%X 99\n'
    r = resolve._parse('star', [text])
    assert (r['ra'], r['dec']) == (10, -20)
    assert r['parallax_mas'] is None
    assert r['ref_epoch'] is None
    assert r['rv_kms'] is None
    assert resolve.propagate(r, 61000) == (10, -20)


def test_unavailable_names_do_not_prevent_coordinate_query(monkeypatch):
    from api.params import parse
    def unavailable(*args, **kwargs):
        raise httpx.ConnectError('offline')
    monkeypatch.setattr(resolve.httpx, 'get', unavailable)
    from api.errors import ApiError
    with pytest.raises(ApiError) as e:
        resolve.resolve('uncached-unknown-source')
    assert e.value.code == 'upstream'
    assert resolve.for_params(parse(ra='10', dec='20'))['center'] == (10, 20)


def test_midpoint_tai_and_utc_are_same_physical_time():
    from Exposure.loader import _read_mjd
    t = Time(60903, format='mjd', scale='utc')
    assert _read_mjd([{'MJD-AVG': t.tai.mjd, 'TIMESYS': 'TAI'}], None) == pytest.approx(t.mjd, abs=1e-10)
    assert _read_mjd([{'MJD-BEG': 1}, {'DATE-AVG': t.isot}], None) == pytest.approx(t.mjd)


def test_fixed_and_tracking_queries_never_share_cache():
    from api import jobs
    from api.params import parse
    p = parse(ra='10', dec='20')
    fixed = {'center': (10, 20), 'motion': None}
    moving = {'center': (10, 20), 'motion': resolve.resolve('barnard')}
    assert jobs._view_key(p, fixed) != jobs._view_key(p, moving)


def test_success_html_is_service_failure_not_unknown_name(monkeypatch):
    from api.errors import ApiError
    monkeypatch.setattr(resolve.httpx, 'get', lambda *args, **kw: httpx.Response(200, text='<html>maintenance</html>'))
    with pytest.raises(ApiError) as error:
        resolve._fetch_sesame('unknown')
    assert error.value.code == 'upstream'


def test_sesame_canonical_identifier_and_missing_motion():
    record = resolve._parse('alias', ['#=Sc=Simbad\n%J 10 20\n%I.0 HD 123\n'])
    assert record['source_id'] == 'HD 123'
    assert record['pmra_masyr'] is None and record['parallax_mas'] is None
    assert resolve.propagate(record, 61000) == (10, 20)


def test_public_and_internal_catalog_aliases():
    for name in ('eps Ind B', 'epsilon_indi_b', 'wise_0855', 'rw_aur', 'r_dra'):
        assert resolve.resolve(name)['coordinate_model'] == 'geocenter_geometric'


def test_sesame_retains_available_rv_and_scientific_notation():
    r = resolve._parse('alias', ['#=Sc=Simbad\n%J 1e1 -2e1\n%P -1e2 +2e2\n%V v -2.05e1 A\n%X 1e-1\n'])
    assert r['rv_kms'] == -20.5
    assert (r['pmra_masyr'], r['pmdec_masyr']) == (-100, 200)
    assert r['errors'] is None and r['ref_epoch'] is None


def test_time_source_scale_and_lower_bound_are_explicit():
    from Exposure.loader import _read_time_info
    t = Time(60903, format='mjd', scale='utc')
    headers = [{'TIMESYS': 'TT'}, {'MJD-AVG': t.tai.mjd, 'TIMESYS': 'TAI'}]
    mjd, source, assumptions = _read_time_info(headers, None)
    assert mjd == pytest.approx(t.mjd, abs=1e-10) and source == 'MJD-AVG' and not assumptions
    assert 'timestamp_not_verified_midpoint' in _read_time_info([{'MJD-BEG': 60903}], None)[2]
    assert 'archive_time_may_be_lower_bound' in _read_time_info([{}], 60903)[2]


@pytest.mark.parametrize('line', ['%J NaN 20', '%J 999 20', '%J 10 1e999'])
def test_malformed_upstream_position_is_not_an_unknown_name(line):
    from api.errors import ApiError
    with pytest.raises(ApiError) as error:
        resolve._parse('real-name', ['#=Sc=Simbad\n'+line+'\n'])
    assert error.value.status == 502 and error.value.code == 'upstream'
