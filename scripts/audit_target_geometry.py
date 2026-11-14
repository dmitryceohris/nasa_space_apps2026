"""Audit API marker axes and geometric margins using saved responses and real FITS headers.
This does not measure registration or source centroids. FITS stay outside the repository.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import sky_motion


def audit(responses, fits_cache):
    native = {}
    for path in fits_cache.glob('cut_*.fits'):
        try:
            header = fits.getheader(path, 'IMAGE')
            native.setdefault((str(header.get('OBSID')), int(header['DETECTOR'])), []).append(header)
        except (OSError, KeyError, ValueError):
            continue
    targets = sky_motion.load_targets()
    results = []
    for path in sorted(responses.glob('[0-7]-*.json')):
        record = json.loads(path.read_text())
        view = record.get('responses', {}).get('view', {})
        obj = record.get('responses', {}).get('resolve', {})
        if not view.get('frames'):
            continue
        grid = WCS(view['grid']['wcs'])
        height, width = view['grid']['height_px'], view['grid']['width_px']
        for frame in view['frames']:
            moment = Time(frame['mjd'], format='mjd', scale='utc')
            coord = (sky_motion.predict_position(targets[obj['name']], moment).geocentric
                     if obj.get('coordinate_model') == 'geocenter_geometric'
                     else SkyCoord(obj['ra'] * u.deg, obj['dec'] * u.deg))
            x, y = [float(v) for v in grid.world_to_pixel(coord)]
            expected = np.array([x, height - 1 - y])
            actual = frame.get('predicted_px')
            error = float(np.max(np.abs(expected - actual))) if actual is not None else None
            check = {'target': record['input'], 'date': frame['date'], 'obs_ids': frame['obs_ids'],
                     'grid_margin_px': min(x, y, width - 1 - x, height - 1 - y),
                     'marker_error_px': error, 'native': []}
            for obs_id in frame['obs_ids']:
                for header in native.get((obs_id, frame['det']), []):
                    wcs = WCS(header)
                    xx, yy = [float(v) for v in wcs.world_to_pixel(coord)]
                    margin = min(xx, yy, header['NAXIS1'] - 1 - xx, header['NAXIS2'] - 1 - yy)
                    check['native'].append({'obs_id': obs_id, 'margin_px': margin,
                                            'background_radius_7px_fits': bool(margin >= 7),
                                            'time_source': next((key for key in ('MJD-AVG','DATE-AVG','MJD-BEG') if key in header), None),
                                            'timesys': header.get('TIMESYS')})
            results.append(check)
    return {'scope':'geometric FITS/WCS audit; no empirical registration or source detection', 'frames':results,
            'marker_axes_pass': bool(results) and all(r['marker_error_px'] is not None and r['marker_error_px'] < .002 for r in results),
            'all_grid_margins_include_7px': bool(results) and all(r['grid_margin_px'] >= 7 for r in results),
            'all_native_margins_include_7px': bool(results) and all(r['native'] and all(n['background_radius_7px_fits'] for n in r['native']) for r in results)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--responses', type=Path, default=Path('data/validation/cold/responses'))
    parser.add_argument('--fits-cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('data/validation/geometry_audit.json'))
    args = parser.parse_args()
    result = audit(args.responses, args.fits_cache)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print({key:value for key,value in result.items() if key!='frames'}, 'frames=',len(result['frames']))
