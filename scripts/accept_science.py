"""Reproducible science audit: saved evidence plus a bounded all-detector archive sample.

Offline by default. --refresh-archive records fresh SIA responses; --per-detector
explicitly bounds downloaded distinct dates, not a claim of archive completeness.
Large FITS stay in cache; JSON retains URLs, SHA-256, rejected measurements and reasons.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
from io import BytesIO
import json
from pathlib import Path
import sys
import warnings

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'src')]

import httpx
import numpy as np
from astropy.io.votable import parse_single_table
from astropy.time import Time

from Exposure import load_exposure
from processing.observations import merge_registry, product_record
from processing.observer import observer_position
from processing.photometry import forced_photometry
from processing.psf import fit_epsf
from processing.targets import load_targets
from spxhack import predict


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def original_selections():
    rows=[]
    for file in ['data/irsa_barnard_20261114/observations.csv','data/physics_observations.csv']:
        with (ROOT/file).open() as stream:
            for r in csv.DictReader(stream):
                rows.append({'product_url':r['access_url'],'obs_id':r['obs_id'],'selection':file})
    return rows


def audit_saved(cache_dir):
    evidence=json.loads((ROOT/'data/validation/science-acceptance/evidence.json').read_text())
    wanted=evidence['frames']+evidence['wise_exact_coordinate_frames']
    hashes={digest(p):p for p in cache_dir.glob('cut_*.fits')}
    rows=[]
    for r in wanted:
        item=dict(r,selection='saved-eight-target-acceptance',path=hashes.get(r['cutout_sha256']))
        rows.append(item)
    return rows


def query_target(t,out,refresh):
    path=out/f'{t.id}-search.vot'
    status={'target':t.id,'source':'saved_current_audit' if path.exists() else 'unavailable'}
    if refresh:
        center=t.at(61359)  # 2026-11-15; broad 0.05 deg covers motion over the survey.
        params=[('COLLECTION','spherex_dr1'),('COLLECTION','spherex_qr3'),
                ('POS',f'CIRCLE {center.ra.deg:.6f} {center.dec.deg:.6f} 0.05'),('DPTYPE','image'),
                ('CALIB','2'),('MAXREC','10000'),('RESPONSEFORMAT','VOTABLE')]
        try:
            response=httpx.get('https://irsa.ipac.caltech.edu/SIA',params=params,timeout=45,follow_redirects=True)
            response.raise_for_status()
            table=parse_single_table(BytesIO(response.content)).to_table(use_names_over_ids=True)
            path.write_bytes(response.content)
            status.update(source='live_sia',url=str(response.url))
        except Exception as error:
            status.update(source='archive_error',error=str(error))
    if not path.exists():
        fallback={'barnard':'data/irsa_barnard_20261114/search_results.vot',
                  'epsindib':'data/physics_evidence/epsilon_indi_b/search_results.vot',
                  'wise0855':'data/physics_evidence/wise_0855/search_results.vot',
                  'rwaur':'data/physics_evidence/rw_aur/search_results.vot',
                  'rdra':'data/physics_evidence/r_dra/search_results.vot'}[t.id]
        saved=ROOT/fallback
        if not saved.exists(): return [],status
        path.write_bytes(saved.read_bytes())
        status.update(source='saved_2026_11_14_after_archive_failure',fallback=fallback)
    table=parse_single_table(path).to_table(use_names_over_ids=True)
    status['response_sha256']=digest(path)
    rows=[]
    for r in table:
        if not str(r['energy_bandpassname']).startswith('SPHEREx-D'): continue
        det=int(str(r['energy_bandpassname']).rsplit('D',1)[1])
        rows.append({'target':t.id,'detector':det,'product_url':str(r['access_url']),
                     'obs_id':str(r['obs_id']),'mjd_hint':float(r['t_min']+r['t_max'])/2,
                     'release':str(r['obs_collection']).removeprefix('spherex_')})
    status.update(n_candidates=len(rows),detectors=dict(Counter(r['detector'] for r in rows)),
                  retrieval_attempt_at_utc=Time.now().utc.isot,scope='metadata coverage search, not exhaustive pixel quality')
    return rows,status


def select_dates(rows,count):
    chosen=[]
    for det in range(1,7):
        sub=sorted([r for r in rows if r['detector']==det],key=lambda r:(r['mjd_hint'],r['release']!='dr1',r['obs_id']))
        # Prefer DR1 for the same exposure; preserve release variants in the registry.
        unique={}
        for r in sub:
            key=r['obs_id']
            if key not in unique or r['release']=='dr1': unique[key]=r
        dates={}
        for r in sorted(unique.values(),key=lambda r:r['mjd_hint']):
            dates.setdefault(int(r['mjd_hint']),r)
        ordered=list(dates.values())
        if count and ordered:
            chosen.extend(ordered[j] for j in np.unique(np.linspace(0,len(ordered)-1,min(count,len(ordered))).astype(int)))
    return chosen


def download(item,cache_dir,targets,allow_download=False):
    t=targets[item['target']]
    c=predict(t,item['mjd_hint'])
    url=item['product_url']+f'?center={c.ra.deg:.12f},{c.dec.deg:.12f}&size=0.05deg&gzip=false'
    path=cache_dir/('science_'+hashlib.sha256(url.encode()).hexdigest()[:24]+'.fits')
    try:
        if not path.exists():
            if not allow_download:
                return dict(item,path=None,cutout_url=url,selection=item.get('selection','six-detector-sample'),
                            error='not_cached_offline')
            response=httpx.get(url,timeout=45,follow_redirects=True)
            response.raise_for_status()
            # Parse before publishing a successful download.
            from astropy.io import fits
            with fits.open(BytesIO(response.content)) as hd:
                if 'IMAGE' not in hd or 'VARIANCE' not in hd: raise ValueError('invalid FITS')
            path.write_bytes(response.content)
        return dict(item,path=path,cutout_url=url,selection=item.get('selection','six-detector-sample'))
    except Exception as error:
        return dict(item,path=None,cutout_url=url,selection='six-detector-sample',error=str(error))


def measure(item,targets,providers=None):
    record={k:v for k,v in item.items() if k not in ('path','wavelength_samples','local_flags_5x5')}
    record['selection']=item['selection']
    if item['path'] is None:
        record.update(status='unavailable',reason=item.get('error','saved_fits_missing'))
        return record
    e=load_exposure(str(item['path']),product_url=item['product_url'])
    record.update(sha256=e.provenance['file_sha256'],mjd=e.mjd,time_source=e.time_source,
                  provenance=e.provenance)
    t=targets.get(item['target'])
    # For the three compatibility stars use the verified saved coordinate, without inventing motion.
    coord=predict(t,e.mjd) if t else e.pixel_to_sky(*item['native_xy'])
    if not e.contains(coord):
        record.update(status='rejected',reason='target_outside_image')
        return record
    p=fit_epsf(e,coord) if e.calibration is not None and e.epsf_data is not None else forced_photometry(e,coord)
    observer=observer_position(e)
    if t and observer['orbit_bound_au'] is not None:
        observer['max_orbital_parallax_arcsec']=observer['orbit_bound_au']*t.parallax/1000
    record.update(photometry=p,observer=observer,status=p['validation']['status'] if p else 'rejected',
                  astrometry={'status':'provisional','reasons':['independent_target_identity_unverified','reference_registration_required']})
    if t:
        from spxhack.native import native_track
        track=native_track([e],t,reference_catalog=(providers or {}).get(t.id))
        native=track.rows[0]
        record['astrometry']=native['astrometry_validation']
        record['native_centroid']={k:native.get(k) for k in ('ra','dec','xy','covariance_arcsec2','association','registration')}
        record['n_independent_detections']=len(native['candidates'])
        record['measurement_scope']='unresolved_system' if t.id in ('epsindib','rwaur') else 'catalogue_position_forced_photometry'
    return record


def run(cache_dir,out,refresh=False,per_detector=0,download_missing=False):
    out.mkdir(parents=True,exist_ok=True);cache_dir.mkdir(parents=True,exist_ok=True)
    targets=load_targets()
    providers={}
    from processing.references import reference_provider
    for t in targets.values():
        path=cache_dir/'references'/f'{t.id}.csv'
        if path.exists(): providers[t.id]=reference_provider(path,t)
    rows=audit_saved(cache_dir)
    registry_rows=original_selections()
    searches=[];selected=[]
    for t in targets.values():
        candidates,status=query_target(t,out,refresh)
        searches.append(status)
        registry_rows.extend(dict(r,selection='all-detector-archive-metadata') for r in candidates)
        selected.extend(select_dates(candidates,per_detector))
        print(t.id,status['source'],len(candidates),'metadata candidates',flush=True)
    if per_detector:
        names={t.name:t.id for t in targets.values()}
        for filename in ('data/irsa_barnard_20261114/observations.csv','data/physics_observations.csv'):
            with (ROOT/filename).open() as stream:
                for r in csv.DictReader(stream):
                    rec=product_record(r['access_url'],obs_id=r['obs_id'])
                    rec.update(target=names.get(r.get('target'), 'barnard'),
                               mjd_hint=float(r.get('mjd') or r['time']),selection=filename)
                    selected.append(rec)
    if selected:
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures=[pool.submit(download,r,cache_dir,targets,download_missing) for r in selected]
            for i,f in enumerate(as_completed(futures)):
                rows.append(f.result())
                if (i+1)%12==0: print('Downloaded',i+1,'/',len(selected),flush=True)
    measured=[]
    for i,item in enumerate(rows):
        try: r=measure(item,targets,providers)
        except Exception as error: r={k:v for k,v in item.items() if k!='path'}|{'status':'processing_error','error':str(error)}
        measured.append(r)
        registry_rows.append({k:r[k] for k in ('product_url','obs_id','detector','selection','sha256','cutout_url','mjd','time_source') if k in r})
        if (i+1)%12==0: print('Measured',i+1,'/',len(rows),flush=True)
    summary={}
    for tid in sorted({r['target'] for r in measured}):
        seq=[r for r in measured if r['target']==tid]
        reasons=Counter(reason for r in seq for reason in (r.get('photometry') or {}).get('quality_reasons',[]))
        summary[tid]={'n_products':len({r['product_url'] for r in seq}),'n_cutouts':len(seq),
                      'photometry_status':dict(Counter(r['status'] for r in seq)),
                      'reasons':dict(reasons),'validated_motion':False,'validated_lightcurve':False,
                      'conclusion':'No validated motion or light curve from this bounded dataset; not an archive-wide impossibility claim'}
    report={'checked_at_utc':Time.now().utc.isot,'sample_per_detector':per_detector,
            'searches':searches,'summary':summary,'measurements':measured,
            'limitations':['finite per-detector date sample','full spectral response unavailable',
                           'spacecraft axes not established by GEOCENTER origin label',
                           'independent reference registration and identity required for motion']}
    reference_status=cache_dir/'references/status.json'
    if reference_status.exists(): report['reference_catalogue_downloads']=json.loads(reference_status.read_text())
    (out/'measurements.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    (out/'observations.json').write_text(json.dumps({'products':merge_registry(registry_rows)},indent=2)+'\n')
    (out/'calibrations.json').write_text((cache_dir/'calibrations/manifest.json').read_text())
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fits-cache',type=Path,default=ROOT/'cache')
    parser.add_argument('--output',type=Path,default=ROOT/'data/validation/science-remediation')
    parser.add_argument('--refresh-archive',action='store_true')
    parser.add_argument('--per-detector',type=int,default=0,choices=range(0,7))
    parser.add_argument('--download',action='store_true',help='Allow missing cutout downloads; otherwise strictly offline for FITS')
    args=parser.parse_args()
    warnings.filterwarnings('ignore')
    run(args.fits_cache,args.output,args.refresh_archive,args.per_detector,args.download)
