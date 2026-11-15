"""Science regressions and independent injection oracles; no archive access."""
import copy

import numpy as np
import pytest
from astropy.io import fits
from scipy.special import erf

from Exposure import exposure_from_hdul
from Exposure.calibration import Calibration, parent_pixel
from processing.astrometry import apply_registration, fit_registration, match_reference_stars
from processing.observations import merge_registry
from processing.observer import observer_position, parallax_factors
from processing.photometry import forced_photometry, light_curve
from processing.psf import fit_epsf
from processing.spectral import compare_bandpasses, compatible_groups
from processing.statistics import constant_flux
from test_exposure import _fake_hdul
from test_processing import CENTER, fake_exposure


def calibrated_exposure(seed=0, x0=20.3, y0=19.8, flux=20000, noise=.01):
    e=fake_exposure(None,61000,noise=noise,seed=seed)
    e.header=fits.Header({"CRPIX1A":1.,"CRPIX2A":1.,"RADESYS":"ICRS"})
    e.time_source="MJD-AVG"
    e.calibration=Calibration(np.full(e.shape,2.),np.full(e.shape,.05),3,np.full(e.shape,6.15**2))
    grid=np.arange(65)/5-6.4
    integ=lambda z:.5*(erf((z+.5)/np.sqrt(2))-erf((z-.5)/np.sqrt(2)))
    model=integ(grid)[:,None]*integ(grid)[None,:]
    table=np.zeros(1,dtype=[('XCENTER','f8'),('YCENTER','f8'),('NSTAR','i4'),('EPSF','f8',(65,65))])
    table[0]=(20,20,100,model)
    e.epsf_data=table; e.epsf_header=fits.Header({'OVSMPX':5,'OVSMPY':5})
    yy,xx=np.indices(e.shape)
    truth=integ(yy-y0)*integ(xx-x0)
    k=e.mjysr_to_ujy_per_pixel()
    e.image=(1+.001*(xx-20)+flux*truth/k+np.random.default_rng(seed).normal(0,noise,e.shape)).astype('f4')
    return e,e.pixel_to_sky(x0,y0)


@pytest.mark.parametrize('flag',['PHANMISS','DICHROIC','SATELLITE_HALO'])
def test_required_flags_mask_pixels(flag):
    e=fake_exposure(CENTER,61000)
    e.flags[20,20]=1<<e.flag_bits[flag]
    assert e.bad_mask(True)[20,20]
    p=forced_photometry(e,CENTER)
    assert not p['ok'] and p['flux_ujy'] is None and p['validation']['status']=='rejected'


def test_contamination_never_contributes_and_missing_core_not_rescaled():
    e=fake_exposure(CENTER,61000)
    e.flags[20,20]=1<<e.flag_bits['GHOST']; e.image[20,20]=1e9
    p=forced_photometry(e,CENTER)
    assert p['flux_ujy'] is None and p['err_ujy'] is None
    assert light_curve([p],1,3)==[]
    assert light_curve([p],1,3,include_diagnostics=True)==[p]


def test_missing_flags_and_bad_shapes_are_explicit():
    h=_fake_hdul(); del h['FLAGS']
    with pytest.warns(UserWarning): e=exposure_from_hdul(h)
    assert not e.flags_available
    p=forced_photometry(e,e.pixel_to_sky(10,10),r_in=3,r_out=6)
    assert not p['ok'] and 'flags_missing' in p['quality_reasons']
    h=_fake_hdul(); h['VARIANCE'].data=np.ones((2,2))
    with pytest.raises(ValueError,match='shape mismatch'): exposure_from_hdul(h)


def test_source_and_fullsample_not_masked_and_source_excluded_from_background():
    e=fake_exposure(CENTER,61000)
    e.flags[20,20]=(1<<e.flag_bits['SOURCE'])|(1<<e.flag_bits['FULLSAMPLE'])
    assert not e.bad_mask(True)[20,20]
    e.flags[20,25]=1<<e.flag_bits['SOURCE']; e.image[20,25]=1e6
    p=forced_photometry(e,CENTER)
    e.image[20,25]=0
    q=forced_photometry(e,CENTER)
    assert p['flux_ujy']==q['flux_ujy']


def test_parent_mapping_calibration_and_detector_guard():
    e,_=calibrated_exposure()
    yy,xx=np.indices((100,100))
    c=Calibration(1+xx*.001+yy*.002,np.full((100,100),.05),3)
    e.header['CRPIX1A']=-9.; e.header['CRPIX2A']=-19.
    e.calibration=c
    assert parent_pixel(e,3,4)==(13,24)
    assert e.wavelength_at(3,4)[0]==pytest.approx(1.061)
    c.detector=2
    with pytest.raises(ValueError,match='detector'): e.wavelength_at(3,4)


def test_product_identity_preserves_releases_versions_and_selection_membership():
    url='https://example.org/ibe/data/spherex/dr1/level2/2025W35_1A/l2b-v27-2026-230/1/level2_2025W35_1A_0432_4D1_spx_l2b-v27-2026-230.fits'
    rows=[{'product_url':url+'?center=1,2','selection':'initial'},
          {'product_url':url,'selection':'physics'},
          {'product_url':url.replace('/dr1/','/qr3/'),'selection':'physics'}]
    products=merge_registry(rows)
    assert len(products)==2
    assert products[0]['selections']==['initial','physics']
    assert products[0]['product_id']!=products[1]['product_id']


@pytest.mark.parametrize('phase',[(20.,20.),(20.3,19.8),(20.49,20.49)])
def test_epsf_recovers_flux_and_subpixel_position(phase):
    e,c=calibrated_exposure(x0=phase[0],y0=phase[1])
    p=fit_epsf(e,c)
    assert p['ok'],p
    assert p['flux_ujy']==pytest.approx(20000,rel=.005)
    assert p['stat_err_ujy']>0
    assert p['err_ujy']**2==pytest.approx(p['stat_err_ujy']**2+(p['flux_ujy']*.03)**2)
    # Start away from the true center, as an independent detector would.
    q=fit_epsf(e,e.pixel_to_sky(20.1,20.1),fit_position=True)
    assert q['ok'],q
    assert q['xy']==pytest.approx(phase,abs=.03)
    assert np.min(np.linalg.eigvalsh(q['covariance_xy_px2']))>0


def test_epsf_masked_core_saturation_and_low_snr_rejected():
    e,c=calibrated_exposure()
    e.flags[20,20]=1<<e.flag_bits['OVERFLOW']
    p=fit_epsf(e,c)
    assert not p['ok'] and 'masked_psf_core' in p['quality_reasons']
    e,c=calibrated_exposure(flux=0,noise=.1)
    p=fit_epsf(e,c)
    assert not p['ok'] and 'low_snr' in p['quality_reasons']


def response_point(wave=2.,mjd=61000):
    return {'obs_id':str(mjd),'wave_um':wave,'bw_um':.05,'mjd':mjd,'flux_ujy':100.,'err_ujy':4.,
            'stat_err_ujy':2.,'systematic_err_ujy':3.,'ok':True,'validation':{'status':'validated'},
            'response':{'wave_um':[wave-.025,wave,wave+.025],'throughput':[0,1,0],'source':'synthetic-oracle'}}


def test_lightcurve_requires_full_responses_and_rejects_ok_false():
    a=response_point(); b=response_point(mjd=61100)
    assert len(light_curve([a,b],1,3))==2
    b['ok']=False
    assert light_curve([a,b],1,3)==[a]
    a.pop('response')
    assert light_curve([a],1,3)==[]
    assert compare_bandpasses(a,a)['status']=='provisional'
    assert compare_bandpasses(response_point(),response_point(2.2))['status']=='rejected'


def test_response_comparison_not_just_center_and_groups_not_transitive():
    a=response_point(); b=copy.deepcopy(a)
    b['response']['throughput']=[1,0,1]
    assert compare_bandpasses(a,b)['status']=='rejected'
    a,b,c=response_point(),response_point(2.0002,61100),response_point(2.0004,61200)
    assert all(compare_bandpasses(x,y)['status']=='validated' for g in compatible_groups([a,b,c]) for x in g for y in g)


def test_constant_flux_uses_correlated_calibration_covariance():
    a,b=response_point(),response_point(mjd=61100); b['flux_ujy']=105
    assert constant_flux([a,b])['status']=='provisional'
    independent=constant_flux([a,b],calibration_correlation=0)
    common=constant_flux([a,b],calibration_correlation=1)
    assert common['chi2']>independent['chi2']
    with pytest.raises(ValueError): constant_flux([a,b],covariance_ujy2=[[1,2],[2,1]])


def test_registration_affine_outlier_covariance_and_insufficient_stars():
    rng=np.random.default_rng(5)
    xy=rng.uniform(-30,30,(30,2))
    ref=xy@np.array([[1.001,.002],[-.001,.999]]).T+[.3,-.5]+rng.normal(0,.01,xy.shape)
    ref[0]+=[8,9]
    cov=np.tile(np.eye(2)*.01**2,(len(xy),1,1))
    r=fit_registration(xy,ref,cov,references_propagated=True)
    assert r['validation']['status']=='validated'
    position,C=apply_registration(r,[0,0],np.eye(2)*.01)
    assert position==pytest.approx([.3,-.5],abs=.03)
    assert np.min(np.linalg.eigvalsh(C))>.009
    assert r['n_references']==29
    assert fit_registration(xy[:5],ref[:5],cov[:5])['validation']['status']=='rejected'
    assert not match_reference_stars([[0,0]],[[.1,0],[.2,0]])


def test_spacecraft_origin_is_not_an_axis_frame_and_parallax_bound():
    e,_=calibrated_exposure()
    e.header.update({'X_SC':7000.,'Y_SC':0.,'Z_SC':0.,'XYZ_SC_SYSTEM':'GEOCENTER'})
    fallback=observer_position(e)
    assert fallback['status']=='provisional'
    precise=observer_position(e,spacecraft_frame='ICRS_ALIGNED')
    assert precise['status']=='validated'
    a=parallax_factors(270,10,fallback['position_au'])
    b=parallax_factors(270,10,precise['position_au'])
    assert np.linalg.norm(a-b)<=fallback['orbit_bound_au']*1.0001


def test_epsf_statistical_errors_have_injection_coverage():
    pulls=[]
    for seed in range(40):
        e,c=calibrated_exposure(seed=seed,noise=.04)
        p=fit_epsf(e,c)
        assert p['ok'],p['quality_reasons']
        pulls.append((p['flux_ujy']-20000)/p['stat_err_ujy'])
    assert abs(np.mean(pulls))<.5
    assert .65<np.std(pulls,ddof=1)<1.4


def test_motion_full_covariance_and_observer_inputs():
    from astropy import units as u
    from astropy.time import Time
    from spxhack import Track, Target, fit_motion
    t=Target('test','injected',150.,20.,0.,0.,0.,0.,Time(60000,format='mjd'))
    times=np.arange(6)*100+60000.
    rows=[]
    for i,mjd in enumerate(times):
        coord=CENTER.spherical_offsets_by(2*(mjd-60000)/365.25*u.arcsec,-(mjd-60000)/365.25*u.arcsec)
        rows.append({'index':i,'mjd':mjd,'xy':[20,20],'ra':float(coord.ra.deg),'dec':float(coord.dec.deg),
                     'err_arcsec':.1,'quality':'good','flags':[],'association':'unique_native',
                     'astrometry_validation':{'status':'validated'}})
    tr=Track(rows,t,6.15)
    C=np.eye(12)*.01+np.ones((12,12))*.002
    obs=np.tile([1.,0.,0.],(6,1))
    motion=fit_motion(tr,parallax='none',covariance_arcsec2=C,observer_barycentric_au=obs)
    assert motion.pm==pytest.approx([2.,-1.],abs=.001)
    assert motion.validation['status']=='validated'
    assert motion.covariance.shape==(4,4)
    assert fit_motion(tr,parallax='none').validation['status']=='provisional'
    with pytest.raises(np.linalg.LinAlgError):
        fit_motion(tr,parallax='none',covariance_arcsec2=-np.eye(12))
    with pytest.raises(ValueError,match='degenerate'):
        fit_motion(tr,parallax='fit',observer_barycentric_au=np.zeros((6,3)))


def test_native_detection_keeps_all_reference_stars_not_only_display_choices():
    from processing.astrometry import native_detections
    e,_=calibrated_exposure(flux=0)
    yy,xx=np.indices(e.shape)
    integ=lambda z:.5*(erf((z+.5)/np.sqrt(2))-erf((z-.5)/np.sqrt(2)))
    for x in (8,20,32):
        for y in (8,20,32):
            e.image+=20000*integ(yy-y)*integ(xx-x)/e.mjysr_to_ujy_per_pixel()
    points=native_detections(e)
    assert len(points)==9
    assert all(p['ok'] and not p['position_assisted'] for p in points)


def test_simultaneous_neighbour_fit_recovers_target_and_exposes_blend():
    e,c=calibrated_exposure()
    yy,xx=np.indices(e.shape)
    integ=lambda z:.5*(erf((z+.5)/np.sqrt(2))-erf((z-.5)/np.sqrt(2)))
    e.image+=15000*integ(yy-19.8)*integ(xx-23.5)/e.mjysr_to_ujy_per_pixel()
    bad=fit_epsf(e,c)
    assert not bad['ok']
    good=fit_epsf(e,c,neighbours=[(23.5,19.8)])
    assert good['ok'],good
    assert good['flux_ujy']==pytest.approx(20000,rel=.01)
    assert good['neighbour_count']==1


def test_independent_native_track_registers_against_eight_reference_stars():
    from astropy.time import Time
    from spxhack import Target, native_track
    e,_=calibrated_exposure(flux=0)
    yy,xx=np.indices(e.shape)
    integ=lambda z:.5*(erf((z+.5)/np.sqrt(2))-erf((z-.5)/np.sqrt(2)))
    refs=[]
    for x in (8,20,32):
        for y in (8,20,32):
            e.image+=20000*integ(yy-y)*integ(xx-x)/e.mjysr_to_ujy_per_pixel()
            if (x,y)!=(20,20): refs.append(e.pixel_to_sky(x,y))
    e.header.update({'X_SC':7000.,'Y_SC':0.,'Z_SC':0.,'XYZ_SC_SYSTEM':'GEOCENTER','XYZ_SC_FRAME':'ICRS_ALIGNED'})
    t=Target('synthetic','test',150.,20.,0.,0.,0.,0.,Time(60000,format='mjd'))
    provider=lambda exp:{'sky':refs,'covariance_px2':np.tile(np.eye(2)*1e-8,(8,1,1)),'source':'independent synthetic truth'}
    row=native_track([e],t,reference_catalog=provider)[0]
    assert row['astrometry_validation']['status']=='validated',row
    assert row['xy']==pytest.approx([20,20],abs=.03)
    assert row['registration']['n_references']==8
    assert row['association']=='unique_native'
