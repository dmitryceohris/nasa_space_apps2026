# Running and deploying SPHEREx Viewer

The main page and API are deployed on Vercel. Completed results and images are stored in PostgreSQL on Supabase. [Open SPHEREx Sky Viewer](https://spherex-sky-viewer-zeta.vercel.app/). Deployment settings and coordinate links are described in “Vercel and Supabase”.

## Local operation

```sh
uv sync --locked
PREWARM=0 bash start.sh
```

Open `http://127.0.0.1:8000/`. To select a port and external interface:

```sh
HOST=0.0.0.0 PORT=8000 PREWARM=0 bash start.sh
```

`CENTER_MJD` is optional. It fixes the query-centre date for a reproducible demo; in normal operation the centre for a name is calculated for the current date. Entered RA/Dec coordinates are always fixed.

## Container

The `Dockerfile` uses Python 3.12.13, uv 0.11.29, and `uv.lock`. The application runs as a separate user; `/cache` is intended for a persistent volume. Personal environments, Git, runtime cache, FITS files, and validation responses are excluded.

```sh
docker compose up --build -d
```

Open `http://127.0.0.1:8000/`. The named `spherex-cache` volume retains data across application restarts. `docker compose down` leaves the volume intact. On the chosen host, mount persistent storage at `/cache`, provide `PORT` if necessary, and configure an external address and HTTPS in front of the application.

One application process is used. The queue and locks are in memory: multiple processes sharing a cache directory require a separate shared queue and locks.

## Data and load limits

By default, images are searched in the live archive. Set `SPHEREX_SOURCE=tables` for an intentionally limited set from local tables; this mode uses separate cache keys. In `archive` mode, a local table does not hide new observations.

| Variable | Default | Meaning |
| --- | ---: | --- |
| `CACHE_DIR` | `cache/`, `/cache` in container | FITS, metadata, PNG and results |
| `PREWARM` | 0 in container | Do not automatically start the five targets |
| `RESULT_CACHE_TTL` | 3600 seconds | Lifetime of a completed result |
| `PARTIAL_CACHE_TTL` | 60 seconds | Retry missed downloads after `partial`, `download_failed`, or `unusable_data` |
| `EMPTY_CACHE_TTL` | 60 seconds | Refresh interval for an empty response |
| `EXPOSURE_CACHE_TTL` | 3600 seconds | Refresh interval for the exposure list |
| `TAP_PROBE_TIMEOUT` | 12 seconds | One TAP attempt before automatic fallback to official IRSA SIA2 search |
| `TAP_RETRY_SECONDS` | 300 seconds | Delay before retrying TAP after failure; SIA2 continues serving requests |
| `MAX_PENDING_JOBS` | 16 | Queued and running jobs (saved errors excluded); subsequent requests get HTTP 503 `busy` |
| `MAX_JOB_SECONDS` | 1800 seconds | Total wait/run limit, checked between stages and downloads |
| `JOB_RETENTION_SECONDS` | 300 seconds | Retention of forgotten completed errors |

Four jobs run concurrently. On timeout, remaining queued downloads are canceled; HTTP calls already started finish under their own network limits. While they finish, a repeated request does not create a second job with the same key. IRSA API calls use three 45-second attempts with a bounded delay. The spectral library retains its own download limit and retries. Damaged FITS and incomplete JSON are not published as completed results.

In the Vercel entry point (`app.py`), `SPHEREX_EPHEMERAL_STORAGE=1` enables one heavy calculation per instance, a separate working directory, per-epoch downloads, and cleanup after success, error, or timeout. PNGs are removed after the Supabase transaction is confirmed. Working directories from interrupted calls are removed on the next job. Inactive temporary-cache files are evicted by last use; open files are protected.

All Vercel application directories are isolated inside `/tmp/spherex-runtime`. Platform-provided XDG values are not used as cleanup roots: Vercel service files and installed dependencies are not affected. `SPHEREX_TMP_MAX_BYTES=350000000` limits the application's total temporary cache, including calibrations and library caches. `SPHEREX_TMP_MIN_FREE_BYTES=75000000` reserves free filesystem space. If space cannot be safely freed, the API returns `503`, `error.code=storage_full`, and `Retry-After: 5`, rather than a completed result without pixels. The `scratch job finished` log contains peak managed-file usage and remaining free space. On Vercel, calibrations are checked by SHA-256 and held in a bounded memory cache: the platform uses a substantial part of `/tmp` for its own dependencies, so an additional 50 MB of calibrations per detector cannot be stored there.

A local `api.main:app` run keeps the normal FITS cache and does not enable archive cleanup. Vercel requires a configured `SPHEREX_DATABASE_URL`; completed frames and GIFs are read from the database independently of temporary-directory lifetimes.

## Verification after deployment

Replace `BASE` with the external address:

```sh
.venv/bin/python scripts/validate_targets.py --base BASE --images --timeout 300 --equivalent-coordinates --output data/validation/published
```

Then, in a normal browser on another device:

1. Open the empty home page; the form should not initiate image downloads.
2. Check the five presets and the exact WISE pair `08h55m10.831s, -7d14m42.53s`, field 6′, D6.
3. Switch dates and check PNG availability and warnings.
4. Open manifest links for three additional stars and equivalent coordinates.
5. Check a repeated query, an unknown name, an invalid or incomplete coordinate pair, and archive failure.
6. Check `/dev` and `/dev/` if the test page remains available.

This list and the container alone do not establish that deployment has occurred. For real-frame checks and quality limitations see [TARGET_VALIDATION.md](TARGET_VALIDATION.md).

## API acceptance follow-up

Five fixes by Asanali from `b55b05e` were carried over onto current main, including invalid-spectrum cache refresh and excluding saved errors from the queue limit. Full report: [reviews/api-acceptance.md](reviews/api-acceptance.md).

Internal exceptions were additionally hidden from public responses and background jobs; `Retry-After` was added for 503 responses and expensive API requests were limited with 429. `API_RATE_LIMIT_PER_MINUTE=120` is one client's combined budget across all expensive endpoints; normal polling every 2–2.5 seconds fits within it. When deploying behind a proxy, configure trusted proxy addresses in Uvicorn (`FORWARDED_ALLOW_IPS`) based on the actual deployment network. Without this, the proxy may be treated as one client. Do not trust forwarded headers from arbitrary external addresses.

A single local application process and an unverified live cold spectral run remain limitations. This patch does not select hosting or validate scientific interpretation.

## Vercel and Supabase

The public viewer accepts ICRS coordinates in degrees or hours, minutes and seconds. Example archive query: `/?ra=133.79512917&dec=-7.24514722&det=6&fov=6&mode=single&max_per_epoch=1`. The same URL works as an external iframe's `src`. **Share or embed this view** provides a link and HTML with the current parameters. The field of view is specified in arcminutes.

By default the viewer queries IRSA. PostgreSQL results and exposure lists are refreshed on the next request after their TTL expires. The separate `snapshot=2026-11-15` parameter opens saved frames for the five presets or the exact WISE coordinate selection. These support `single`, `coadd`, `diff`, and `trail`, with `max_per_epoch=6` (default) or `1`. Results were recalculated from verified original FITS, and PNGs and responses are included in the deployment: viewing requires no database, new IRSA search, or background processing. Manifest and parameters: `frontend/challenge-demo/snapshots/manifest.json`; reproduce with `python scripts/export_viewer_snapshots.py` using the original local FITS cache. All network access during export is prohibited. Scientific interpretation limitations and quality masks are retained.

When a user selects one of the five presets, the interface selects saved observations if the user has not explicitly selected a source. `Latest IRSA` remains available for a new search and arbitrary parameters. Saved frames load from the current deployment even when an external server is specified with `?api=`; spectrum and light curve use the archive API.

The slider displays the selected field of view in arcminutes. After resizing or changing detector, the viewer automatically selects `Latest IRSA` and loads with the new parameters; processing may take time. A later preset selection preserves that source. When manually returning to saved observations, field size and detector are restored from the manifest for the selected object or exact WISE coordinate selection.

For the next release, set `SPHEREX_COLLECTIONS` in Vercel (comma-separated SIA2 names in preference order) and redeploy. Current values are `spherex_dr1,spherex_qr3`. This setting applies to fallback SIA2 search and spectra. The primary TAP search uses the current IRSA table. Archive-cache keys include the collection list. When changing calibrations, update the verified manifest and compatible pipeline versions. [IRSA collections and release cadence](https://irsa.ipac.caltech.edu/spherex-archive-documentation/spherex-data-access/).

The server uses the Vercel secret `SPHEREX_DATABASE_URL` for the `spherex` schema. For a new database, apply the migration in `supabase/migrations/`, grant the application access to that schema, and run `python scripts/seed_supabase_demo.py` with the same connection. Frames and masks are stored in PostgreSQL together with results.
