"""Native-image astrometry; no target-catalogue position is used to pick a candidate."""
import numpy as np

from processing.astrometry import apply_registration, native_detections, register_catalogue
from processing.observer import observer_position
from processing.validation import quality


def native_track(exposures, where, *, reference_catalog=None):
    """Independent detections + free ePSF centroids + optional propagated reference catalogue.

    reference_catalog(exp) -> {'sky': per-epoch SkyCoord sequence excluding target,
      'covariance_px2': (N,2,2), 'source': URL or local manifest identifier}.
    The supplied positions must include proper motion and relevant parallax.
    No catalogue means provisional centroids, never verified astrometry.
    """
    from .data import target
    from .measure import Track, predict
    t=target(where)
    if reference_catalog is None:
        from pathlib import Path
        path=Path(__file__).resolve().parents[2]/'cache/references'/f'{t.id}.csv'
        if path.exists():
            from processing.references import reference_provider
            reference_catalog=reference_provider(path,t)
    rows=[]
    for i,exp in enumerate(sorted(exposures,key=lambda e:e.mjd)):
        detections=native_detections(exp)
        refs=reference_catalog(exp) if reference_catalog is not None else None
        reg=register_catalogue(exp,detections,refs['sky'],refs.get('covariance_px2')) if refs else None
        # Remove known reference stars, without using the target's predicted position.
        candidates=[p for p in detections if p.get('ok')]
        if refs:
            candidates=[p for p in candidates if all(
                np.hypot(*(np.asarray(p['xy'])-exp.sky_to_pixel(c)))>1.5 for c in refs['sky'])]
        p=candidates[0] if len(candidates)==1 else None
        row={'index':i,'mjd':exp.mjd,'obs_id':exp.obs_id,'date':str(exp.mjd),
             'xy':None,'x':None,'y':None,'ra':None,'dec':None,'err_x':None,'err_y':None,
             'err_arcsec':None,'quality':None,'flags':[],'association':'unique_native' if p else None,
             'status':'ok' if p else 'ambiguous' if candidates else 'none',
             'reason':None if p else 'no_unique_independent_detection',
             'method':'native_epsf_fit','candidates':detections,'provenance':exp.provenance,
             'registration':reg,'observer':observer_position(exp),'offset_px':None,
             'pred_xy':list(exp.sky_to_pixel(predict(t,exp.mjd)))}
        if p:
            xy,cov=p['xy'],p['covariance_xy_px2']
            if reg and reg['matrix'] is not None:
                xy,cov=apply_registration(reg,xy,cov)
            c=exp.pixel_to_sky(*xy)
            # Numerical Jacobian to local east/north; preserves rotated/distorted WCS.
            columns=[]
            for delta in ([.001,0],[0,.001]):
                shifted=exp.pixel_to_sky(*(np.array(xy)+delta))
                de,dn=c.spherical_offsets_to(shifted)
                columns.append([de.arcsec/.001,dn.arcsec/.001])
            J=np.asarray(columns).T
            skycov=J@np.asarray(cov)@J.T
            row.update(xy=xy,x=xy[0],y=xy[1],ra=float(c.ra.deg),dec=float(c.dec.deg),
                       err_x=float(np.sqrt(cov[0][0])),err_y=float(np.sqrt(cov[1][1])),
                       err_arcsec=float(np.sqrt(np.trace(skycov)/2)),covariance_arcsec2=skycov.tolist(),
                       quality='good',offset_px=float(np.linalg.norm(np.asarray(xy)-row['pred_xy'])))
        row['astrometry_validation']=quality({
            'poor_centroid':p is not None,
            'astrometric_registration_failed':bool(reg and reg['validation']['status']=='validated'),
            'reference_provenance_missing':bool(refs and refs.get('source')),
            'observer_unverified':row['observer']['status']=='validated',
            'unverified_target_identity':p is not None,
            'unresolved_components':t.id not in ('epsindib','rwaur'),
        },['poor_centroid'])
        rows.append(row)
    return Track(rows,t,6.15)
