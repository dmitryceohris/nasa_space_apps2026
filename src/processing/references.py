"""Gaia reference stars: saved source provenance and per-epoch covariance propagation."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import httpx
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from .observer import observer_position, parallax_factors

PARAMETERS=('ra','dec','parallax','pmra','pmdec')
CORRELATIONS=('ra_dec','ra_parallax','ra_pmra','ra_pmdec','dec_parallax','dec_pmra','dec_pmdec',
              'parallax_pmra','parallax_pmdec','pmra_pmdec')


def fetch_references(target,directory,radius_deg=.15):
    """Save exact ESA response; no fabricated replacement if the archive is unavailable."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    columns=['source_id','ra','dec','ref_epoch','parallax','pmra','pmdec','ruwe','phot_g_mean_mag']
    columns += [p+'_error' for p in PARAMETERS]+[p+'_corr' for p in CORRELATIONS]
    c=target.at(61359)
    query=(f"SELECT {', '.join(columns)} FROM gaiadr3.gaia_source WHERE phot_g_mean_mag<18 "
           f"AND ruwe<1.4 AND astrometric_params_solved IN (31,95) AND duplicated_source='false' "
           f"AND 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{c.ra.deg:.7f},{c.dec.deg:.7f},{radius_deg}))")
    response=httpx.get('https://gea.esac.esa.int/tap-server/tap/sync',
                       params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},
                       timeout=45,follow_redirects=True)
    response.raise_for_status()
    rows=list(csv.DictReader(response.text.splitlines()))
    if not rows or not all(c in rows[0] for c in columns): raise ValueError('invalid_or_empty_gaia_reference_response')
    path=directory/f'{target.id}.csv';path.write_bytes(response.content)
    manifest={'url':str(response.url),'sha256':hashlib.sha256(response.content).hexdigest(),
              'frame':'ICRS','reference_epoch_scale':'TCB','radius_deg':radius_deg,
              'n_stars':len(rows),'source':'ESA Gaia DR3','target':target.id}
    path.with_suffix('.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return path


def reference_provider(path,target):
    """Return native_track's catalogue callback using the saved five-parameter covariance."""
    path=Path(path)
    manifest=json.loads(path.with_suffix('.json').read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['sha256']:
        raise ValueError('reference_catalogue_checksum_mismatch')
    rows=list(csv.DictReader(path.open()))
    def provide(exp):
        from astropy.time import Time
        observer=observer_position(exp)
        target_c=target.at(exp.mjd)
        expected=[];pixelcov=[];ids=[]
        for row in rows:
            try:
                a,d,pmra,pmdec,plx,epoch=[float(row[k]) for k in ('ra','dec','pmra','pmdec','parallax','ref_epoch')]
                dt=(Time(exp.mjd,format='mjd',scale='utc').tcb.jyear-epoch)
                c=SkyCoord(a*u.deg,d*u.deg).spherical_offsets_by(pmra*dt*u.mas,pmdec*dt*u.mas)
                factors=parallax_factors(c.ra.deg,c.dec.deg,observer['position_au'])
                c=c.spherical_offsets_by(factors[0]*max(plx,0)*u.mas,factors[1]*max(plx,0)*u.mas)
                # Exclude the moving target and its immediate neighbourhood.
                if c.separation(target_c).arcsec<6*6.15 or not exp.contains(c,margin_px=5): continue
                errors=np.array([float(row[p+'_error']) for p in PARAMETERS])
                if not np.isfinite(errors).all() or np.any(errors<=0): continue
                R=np.eye(5)
                for i in range(5):
                    for j in range(i+1,5): R[i,j]=R[j,i]=float(row[f'{PARAMETERS[i]}_{PARAMETERS[j]}_corr'])
                C=R*np.outer(errors,errors)
                if np.min(np.linalg.eigvalsh(C))<=0: continue
                J=np.array([[1,0,factors[0],dt,0],[0,1,factors[1],0,dt]])
                skycov=J@C@J.T/1e6
                # Bound unverified low-Earth-orbit parallax instead of inventing axis orientation.
                if observer['status']!='validated':
                    if observer['orbit_bound_au'] is None: continue
                    skycov+=np.eye(2)*(observer['orbit_bound_au']*max(plx,0)/1000)**2
                xy=np.asarray(exp.sky_to_pixel(c))
                columns=[]
                for east,north in ((.001,0),(0,.001)):
                    shifted=c.spherical_offsets_by(east*u.arcsec,north*u.arcsec)
                    columns.append((np.asarray(exp.sky_to_pixel(shifted))-xy)/.001)
                P=np.asarray(columns).T
                expected.append(c);pixelcov.append(P@skycov@P.T);ids.append(row['source_id'])
            except (ValueError,KeyError): continue
        return {'sky':expected,'covariance_px2':np.asarray(pixelcov).reshape(-1,2,2),
                'source':manifest['url'],'source_sha256':manifest['sha256'],'source_ids':ids,
                'assumptions':['linear_tangent_propagation','reference_radial_velocity_not_modelled',
                               'nonpositive_reference_parallax_clamped_to_zero']}
    return provide
