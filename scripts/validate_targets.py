"""Small, sequential real API checks. Run against PREWARM=0 server; --images opts into cutouts.
Responses record live-request provenance; the server may use cached archive/FITS data.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import httpx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', default='http://127.0.0.1:8011')
    parser.add_argument('--images', action='store_true')
    parser.add_argument('--timeout', type=int, default=120)
    parser.add_argument('--only', default='barnard,epsindib,wise0855,rwaur,rdra,Proxima Cen,61 Cyg A,Vega')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'data/validation')
    parser.add_argument('--cache-state', default='server cache possible')
    parser.add_argument('--equivalent-coordinates', action='store_true')
    args = parser.parse_args()
    out = args.output
    (out / 'responses').mkdir(parents=True, exist_ok=True)
    rows = []
    for index, name in enumerate(args.only.split(',')):
        result = {'checked_at': datetime.now(timezone.utc).isoformat(), 'input': name,
                  'provenance': 'real HTTP request; ' + args.cache_state, 'responses': {}}
        record = {'target': name, 'status': 'unverified', 'source_id': '', 'epochs': '', 'limitation': '',
                  'checked_at': result['checked_at'], 'fov_arcmin': 6, 'det': 6 if name == 'wise0855' else 3,
                  'ref_epoch': '', 'ref_epoch_format': '', 'ref_epoch_scale': '', 'link': '',
                  'coordinate_check': 'not_requested'}
        try:
            r = httpx.get(args.base + '/api/resolve', params={'target': name}, timeout=args.timeout)
            body = r.json()
            result['responses']['resolve'] = body
            record.update(status='resolved' if r.status_code == 200 else 'resolve_error', source_id=body.get('source_id', ''))
            record.update({key: body.get(key) for key in ('ref_epoch', 'ref_epoch_format', 'ref_epoch_scale')})
            from urllib.parse import urlencode
            record['link'] = '/?' + urlencode({'target': name, 'fov':6, 'det':record['det'], 'mode':'single', 'max_per_epoch':1})
            record['limitation'] = body.get('limitations', str(body.get('error', '')))
            if args.images and r.status_code == 200:
                query = {'target': name, 'mode': 'single', 'max_per_epoch': 1, 'det': 6 if name == 'wise0855' else 3, 'fov': 6}
                deadline = time.monotonic() + args.timeout
                while time.monotonic() < deadline:
                    view = httpx.get(args.base + '/api/view', params=query, timeout=args.timeout)
                    data = view.json()
                    if view.status_code != 202:
                        result['responses']['view'] = data
                        record.update(status=data.get('data_status', data.get('status', 'error')), epochs=len(data.get('frames', [])))
                        break
                    time.sleep(2)
                else:
                    record['status'] = 'image_timeout'
                if args.equivalent_coordinates and body.get('coordinate_model') == 'fixed_catalogue':
                    coords = dict(query)
                    coords.pop('target')
                    coords.update(ra=body['ra'], dec=body['dec'])
                    deadline = time.monotonic() + args.timeout
                    while time.monotonic() < deadline:
                        view = httpx.get(args.base + '/api/view', params=coords, timeout=args.timeout)
                        if view.status_code != 202:
                            fixed = view.json()
                            result['responses']['coordinate_view'] = fixed
                            named = result['responses'].get('view', {})
                            same = [f['obs_ids'] for f in fixed.get('frames', [])] == [f['obs_ids'] for f in named.get('frames', [])]
                            record['coordinate_check'] = 'same_observations' if same and fixed.get('frames') else 'different_or_empty'
                            break
                        time.sleep(2)
                    else:
                        record['coordinate_check'] = 'timeout'
        except Exception as e:
            record.update(status='request_failed', limitation=str(e))
        result['summary'] = record
        (out / 'responses' / f'{index}-{name.lower().replace(" ", "-")}.json').write_text(json.dumps(result, indent=2, ensure_ascii=False))
        rows.append(record)
        print(name, record['status'], flush=True)
    with (out / 'target_cases.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
