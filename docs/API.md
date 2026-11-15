# Scientific validation: November 15, 2026 update

Existing routes are retained. Photometric `points[]` now contain `validation.status` (`validated`, `provisional`, `rejected`), reasons and check results. `ok` means `validated` only; flat graph arrays contain accepted results. Preliminary aperture photometry is not accepted scientific photometry.

`/api/spectrum` adds `quality_counts`, file and calibration provenance, and statistical and systematic ePSF fit errors. `/api/lightcurve` adds `diagnostic_points`, `series` for separate comparable bands, and `spectral_comparison`. A wavelength window without a full response does not establish spectral comparability; science arrays are empty in that case.

`frames[].validation` explicitly separates image viewing from validated astrometry. The result cache version is 12. Detailed rules, Python API changes and real results: [SCIENCE_REMEDIATION.md](SCIENCE_REMEDIATION.md). The original API description is retained below; when scientific quality rules differ, this update takes precedence.

# SPHEREx blink: API

FastAPI backend. Given coordinates or an object name, it returns SPHEREx frames by epoch (PNG and metadata). Data come live from NASA/IPAC IRSA. Positions and proper motions for five presets come from the team table `OBJECT_PARAMETERS.csv` (Gaia DR3 and Kirkpatrick et al. 2021); other names come from SIMBAD via CDS Sesame.

- Local base address: `http://127.0.0.1:8000`
- Install: `uv sync --locked --python 3.12.13` (API dependencies are in `pyproject.toml`; `requirements.txt` remains for pip installation).
- Start: `bash start.sh` (or `.venv/bin/uvicorn api.main:app --port 8000` from the repository root).
- Offline tests: `uv run --offline --locked pytest -q tests`.
- Image processing is not duplicated in `api/`: frames, background, masks, difference, trail and spectrum come from `src/processing` and `src/Exposure`; `api/` handles HTTP, exposure search, 202/200 jobs, cache and target centres.
- CORS is open to all origins.
- Examples below were captured from a running server; complete versions are in `mock/`. The `upstream` and `internal` error examples were composed from code text and not reproduced on the live server.

## Common parameters

The same parameters apply to the page URL and API (`/api/view`, `/api/export.gif`).

| Parameter | Type | Range / values | Default | Example |
|---|---|---|---|---|
| `target` | str | object name (SIMBAD) or preset `id` | none | `Barnard's Star`, `barnard` |
| `ra` | str / float | degrees `0..360` or `08h55m10.83s`, `08:55:10.83` | none | `133.795` |
| `dec` | str / float | degrees `-90..90` or `-07d14m42.5s` | none | `-7.245` |
| `fov` | float | arcminutes, `1..20` | `5` | `3` |
| `det` | int | `1..6` | `3` | `6` |
| `mode` | enum | `single`, `coadd`, `diff`, `trail` | `coadd` | `diff` |
| `max_per_epoch` | int | `1..20` | `6` | `6` |

Supply either `target` or the pair `ra`+`dec`; when both are supplied, `ra`/`dec` take priority. `ra` and `dec` must be supplied together. Axis formats can be mixed, e.g. `ra=133.7951291667` and `dec=-7d14m42.53s`. RA with colons or `h` is in hours, Dec in degrees. RA 360° (24 hours) is normalized to 0°. An explicit pair fixes the direction: a simultaneous `target` name does not enable a motion model.

## Error format

All endpoints use the same error format, including a missing required parameter (`400 bad_param`) and an unknown route (`404 not_found`):

```json
{
  "error": {
    "code": "bad_param",
    "message": "dec должен быть в диапазоне -90..90"
  }
}
```

| HTTP | `code` | When | Example |
|---|---|---|---|
| 400 | `bad_param` | invalid or missing parameter | `mock/error_bad_param.json` |
| 404 | `not_found` | name not found in SIMBAD | `mock/error_not_found.json` |
| 409 | `not_ready` | GIF export before the result is ready | `mock/error_not_ready.json` |
| 502 | `upstream` | IRSA or SIMBAD unavailable after retries | `mock/error_upstream.json` |
| 500 | `internal` | all other cases; traceback is logged | `mock/error_internal.json` |

---

## GET /api/health

Check that the server is alive. No parameters.

```json
{
  "status": "ok",
  "version": "0.1.0"
}
```

## GET /api/targets

List presets for frontend buttons. No parameters. `id` can be passed as `target=<id>`.

```json
[
  {
    "id": "barnard",
    "name": "Barnard's Star",
    "det": 3,
    "fov": 3,
    "note": "самое быстрое собственное движение, ~10.4″/год"
  },
  {
    "id": "epsindib",
    "name": "eps Ind B",
    "det": 3,
    "fov": 5,
    "note": "пара коричневых карликов, ищите метан"
  },
  {
    "id": "wise0855",
    "name": "WISE 0855-0714",
    "det": 6,
    "fov": 3,
    "note": "виден только на 4–5 мкм"
  },
  {
    "id": "rwaur",
    "name": "RW Aur",
    "det": 2,
    "fov": 3,
    "note": "переменная звезда T Tauri"
  },
  {
    "id": "rdra",
    "name": "R Dra",
    "det": 2,
    "fov": 3,
    "note": "мирида, период ~245 дней"
  }
]
```

`WISE 0855-0714` is not a SIMBAD name, so presets (including `wise0855` and `WISE 0855-0714`) are resolved from the team table without querying SIMBAD.

## GET /api/resolve

Resolve a name to coordinates and proper motion.

| Parameter | Required | Description |
|---|---|---|
| `target` | yes | object name or preset `id` |

Presets (by `id` or name) come from `OBJECT_PARAMETERS.csv`: `source` names the catalogue, `ref_epoch` is the coordinate epoch (2016.0 for Gaia DR3), and `rv_kms` is included. Other names go to SIMBAD with coordinates at `ref_epoch` 2000.0; proper motion and parallax may be `null` if SIMBAD does not know them. SIMBAD responses are cached indefinitely (memory and `cache/resolve.json`). Unknown name: `404`.

```json
{
  "name": "Barnard's Star",
  "id": "barnard",
  "ra": 269.44850252543836,
  "dec": 4.739420051112412,
  "pmra_masyr": -801.5509783684709,
  "pmdec_masyr": 10362.394206546573,
  "parallax_mas": 546.975939730948,
  "rv_kms": -110.46822,
  "ref_epoch": 2016.0,
  "source": "Gaia DR3 / gaiadr3.gaia_source"
}
```

## GET /api/view

The main endpoint returns a series of frames by epoch and works as “get or start”:

1. A cached result returns `200` immediately with the full JSON.
2. An already running job returns `202` with progress.
3. Otherwise a background job starts and returns `202`.

The frontend polls the same URL every 2 seconds until `200`. Parameters: see “Common parameters”. The field centre is calculated for today's date (proper motion); `CENTER_MJD` fixes the date.

### 202: job in progress

`stage` passes through `query`, `download`, `process`, `render`.

```json
{
  "status": "running",
  "job_id": "ab856bfa",
  "stage": "query",
  "done": 0,
  "total": 0,
  "message": "Поиск снимков"
}
```

```json
{
  "status": "running",
  "job_id": "ab856bfa",
  "stage": "download",
  "done": 1,
  "total": 18,
  "message": "Скачано 1 из 18 снимков"
}
```

### 200: ready

Example: `/api/view?target=barnard&fov=3&det=3`. Only the first of five frames is shown for brevity (full response in `mock/view_ready.json`).

```json
{
  "status": "ready",
  "query": {
    "ra": 269.44607270990747,
    "dec": 4.770723702300067,
    "fov": 3.0,
    "det": 3,
    "mode": "coadd",
    "target": "Barnard's Star"
  },
  "grid": {
    "width_px": 29,
    "height_px": 29,
    "pixel_scale_arcsec": 6.15,
    "wcs": {
      "CTYPE1": "RA---TAN",
      "CTYPE2": "DEC--TAN",
      "CRVAL1": 269.44607270990747,
      "CRVAL2": 4.770723702300067,
      "CRPIX1": 15.0,
      "CRPIX2": 15.0,
      "CDELT1": -0.0017083333333333334,
      "CDELT2": 0.0017083333333333334
    }
  },
  "frames": [
    {
      "index": 0,
      "label": "Эпоха 1",
      "mjd": 60908.7776,
      "date": "2025-08-21",
      "mjd_span": [
        60902.8654,
        60914.4848
      ],
      "det": 3,
      "wave_um": 2.0725,
      "n_exp": 5,
      "image_url": "/cache/6eee915f7c97b2db/f0.png",
      "flags_url": "/cache/6eee915f7c97b2db/f0_flags.png",
      "predicted_px": [
        13.839,
        16.075
      ],
      "measured_px": [
        14.025,
        15.943
      ],
      "measured_candidates": [
        {
          "kind": "largest+central",
          "xy": [
            14.025,
            15.943
          ],
          "area_px": 146
        }
      ],
      "obs_ids": [
        "2025W33_2C_0321_1",
        "2025W34_1B_0164_3",
        "2025W34_2B_0030_3",
        "2025W35_1A_0008_3",
        "2025W35_1A_0484_2"
      ]
    }
  ],
  "trail_url": null,
  "motion": {
    "pm_reference_arcsec_yr": [
      -0.8016,
      10.3624
    ],
    "source": "Gaia DR3 / gaiadr3.gaia_source",
    "pm_measured_arcsec_yr": null
  },
  "warnings": [
    "Не удалось получить или обработать 1 из 23 снимков",
    "λ в точке цели меняется между эпохами (1.79-2.18 мкм): у объектов с крутым спектром остатки в разностном кадре могут быть не движением",
    "В центре поля есть насыщенные или загрязнённые пиксели (слой флагов)"
  ],
  "credits": "NASA SPHEREx via NASA/IPAC IRSA"
}
```

Response fields:

| Field | Meaning |
|---|---|
| `query` | Actual parameters; `ra`/`dec` are the field centre (object position today) |
| `grid` | Common TAN grid for all frames; `wcs` converts a click to RA/Dec |
| `frames` | Epochs sorted by `mjd`. An epoch groups images without gaps greater than 30 days; continuous series over 20 days are split into windows (`src/processing/epochs.py`), so there may be more frames than half-year passes |
| `frames[].n_exp`, `wave_um` | Number of exposures in the frame; wavelength (µm) at the target on the frame date from spectral WCS (the linear filter causes it to vary across epochs); `null` if unavailable |
| `frames[].image_url`, `flags_url` | Frame PNG (in `diff`, a relative difference residual: red brighter than template, blue darker) and transparent flags layer. Zodiacal light (ZODI extension) is subtracted. No-data pixels (outside coverage, HARD_BAD_FLAGS, and isolated outliers without SOURCE rejected on native pixels before interpolation) are dark purple `(72, 36, 110)`, a colour absent from the grey sky (R=G=B); saturated cores are white. The flags layer marks contamination and saturation (CROSSTALK, GHOST, STREAK, BLOOM, SNOWBALL, HALO, OVERFLOW, NONLINEAR); SOURCE itself is never painted |
| `frames[].predicted_px` | Predicted target position at the frame date (proper motion and parallax); `null` for `ra`/`dec` without a name |
| `frames[].measured_px` | Centroid of one algorithm-selected candidate when `measurement.status == "ok"`, quality is not `poor`, and distance from an available prediction is ≤1 px (6.15″); always `null` for an out-of-frame prediction. Without a motion model it remains a fixed-field candidate. This threshold is a protective heuristic, not a verified registration error or identity. In diff only the positive lobe is used; with a negative lobe within 5 px only candidates remain (`overlaps_negative_residual`) |
| `frames[].measured_candidates` | Candidates with formal noise error and quality diagnostics; that error does not fully describe registration, PSF and blending |
| `frames[].measurement` | Frame measurement result: `status` (`ok` / `ambiguous` / `none`), `reason` (why no single centre), `quality` (`good` / `fair` / `poor`), `flags`, `err_px` (noise error `[x, y]` in pixels), `n_rejected_hot` (isolated hot pixels rejected), `n_low_snr` (spots with S/N < 5, not selected), `target_association` (`no_prediction` / `no_measurement` / `prediction_outside_field` / `consistent_with_prediction` / `far_from_prediction`), `offset_from_prediction_px` |
| `frames[].span_days` | Epoch duration (last minus first image), days; at most `processing.epochs.max_span_days` |
| `frames[].bg_offsets` | Residual background per exposure after ZODI subtraction, MJy/sr (`null` if too few pixels and no subtraction) |
| `frames[].diff_noise_scale` | `diff`/`trail` only: ratio of actual to propagated difference-frame noise (usually < 1 due to reprojection smoothing); the frame is divided by this number |
| `processing` | How frames and centres were obtained: `epochs` (epoch limits), `background`, `difference` (`null` outside `diff`/`trail`), `measurement` (centroid method, parameters, error model, quality grades, `limitations`, `not_implemented`). Shown in “How positions are measured” |
| `frames[].obs_ids` | IRSA observation identifiers |
| `trail_url` | Motion-trail PNG: epoch residuals coloured by time (older purple, newer yellow) over the template. Populated in `diff` and `trail` with at least two epochs, otherwise `null`. Pixels without data in all epochs use the same dark purple as frames |
| `motion` | `pm_reference_arcsec_yr`: SIMBAD proper motion in arcsec/year (`[RA·cos(Dec), Dec]`); `pm_measured_arcsec_yr` remains `null` (measurement from centres not implemented); entirely `null` if motion is unknown |
| `warnings` | Warning strings (e.g. some images did not download) |

All pixel coordinates (`predicted_px`, `measured_px`, `measured_candidates[].xy`) use PNG coordinates: zero-based `[x, y]`, with `y` increasing downward.

### diff and trail modes

- `diff`: `image_url` contains the frame minus the median of the OTHER epochs (excluding the frame itself), normalized by robust residual scale (MAD), `robust_residual_scale`. This is a relative scale, **not calibrated statistical significance in σ**: interpolated-pixel variance, correlations, background and PSF are not fully propagated. A moving object may produce a positive lobe at its current position and a negative lobe at the reference position. Lobes are detected separately (`measured_candidates[].sign`); only the positive one is used for centroid selection. The pair alone does not prove motion. `trail_url` contains a colour visualization of positive residuals.
- `trail`: frames are as in `coadd`, plus `trail_url`.
- Both modes need at least 2 epochs; otherwise `warnings` explains the limitation (for `diff`, no frames are produced).

Example `/api/view?target=wise0855&fov=3&det=6&mode=diff` (first of four frames; full response in `mock/view_diff_trail.json`):

```json
{
  "status": "ready",
  "query": {
    "ra": 133.75778760470473,
    "dec": -7.242024493644443,
    "fov": 3.0,
    "det": 6,
    "mode": "diff",
    "target": "WISE 0855-0714"
  },
  "grid": {
    "width_px": 29,
    "height_px": 29,
    "pixel_scale_arcsec": 6.15,
    "wcs": {
      "CTYPE1": "RA---TAN",
      "CTYPE2": "DEC--TAN",
      "CRVAL1": 133.75778760470473,
      "CRVAL2": -7.242024493644443,
      "CRPIX1": 15.0,
      "CRPIX2": 15.0,
      "CDELT1": -0.0017083333333333334,
      "CDELT2": 0.0017083333333333334
    }
  },
  "frames": [
    {
      "index": 0,
      "label": "Эпоха 1",
      "mjd": 60800.746,
      "date": "2025-05-05",
      "mjd_span": [
        60800.7007,
        60800.8351
      ],
      "det": 6,
      "wave_um": 4.9436,
      "n_exp": 3,
      "image_url": "/cache/f640d60517f96c47/f0.png",
      "flags_url": "/cache/f640d60517f96c47/f0_flags.png",
      "predicted_px": [
        11.985,
        14.167
      ],
      "measured_px": null,
      "measured_candidates": [
        {
          "kind": "largest",
          "xy": [
            20.406,
            10.047
          ],
          "area_px": 25
        },
        {
          "kind": "central",
          "xy": [
            13.511,
            13.567
          ],
          "area_px": 4
        }
      ],
      "obs_ids": [
        "2025W19_1B_0116_1",
        "2025W19_1B_0116_2",
        "2025W19_1B_0138_1"
      ]
    }
  ],
  "trail_url": "/cache/f640d60517f96c47/trail.png",
  "motion": {
    "pm_reference_arcsec_yr": [
      -8.1237,
      0.6732
    ],
    "source": "Kirkpatrick et al. (2021), Table 5 / VizieR J/ApJS/253/7/plx",
    "pm_measured_arcsec_yr": null
  },
  "warnings": [
    "λ в точке цели меняется между эпохами (4.69-4.94 мкм): у объектов с крутым спектром остатки в разностном кадре могут быть не движением",
    "В центре поля есть насыщенные или загрязнённые пиксели (слой флагов)"
  ],
  "credits": "NASA SPHEREx via NASA/IPAC IRSA"
}
```

### Object centre: `measured_candidates`

The algorithm does not know object size or type in advance. It chooses two candidates:
- `largest`: largest object in the frame;
- `central`: object beneath the frame centre, or the nearest to it;
- `largest+central`: the same object fulfils both roles, so there is one candidate.

If there are no candidates, the list is empty. If `largest` and `central` differ or two objects are equally close to the centre, `measurement.status` is `ambiguous` and `measured_px` is `null`; the server does not guess.

Method (`processing.measurement`): robust background and noise (median/MAD with 3σ clipping), σ=1 pixel smoothing, a 3σ threshold on valid pixels only, connected regions, and a flux-weighted first-moment centroid. No-data pixels (NaN, inf) belong to an object only if surrounded by it (saturated core) and receive the nearest valid neighbour's value, so they do not pull the centre; a mask at the frame edge or next to an object does not become an object.

Candidate fields: `kind`, `sign` (`1`/`-1` in `diff`, otherwise `null`), `xy`, `err_px`, `area_px`, `snr`, `quality`, `flags`. Flags: `masked_pixels`, `saturated_core`, `touches_nodata` and `touches_edge` (centre shifted toward the visible part, quality `poor`), `elongated` (touching objects or trail, quality `fair`). Spots with S/N < 5 are not selected: a 3σ threshold finds them in roughly one in three pure-noise frames. PSF fitting, deblending and neural networks are not used (next stage).

Example (eps Ind B, first of four frames; the rest are in `mock/view_ready_two_candidates.json`):

```json
{
  "status": "ready",
  "query": {
    "ra": 331.0984070266234,
    "dec": -56.801257162911156,
    "fov": 5.0,
    "det": 3,
    "mode": "coadd",
    "target": "Epsilon Indi B"
  },
  "grid": {
    "width_px": 49,
    "height_px": 49,
    "pixel_scale_arcsec": 6.15,
    "wcs": {
      "CTYPE1": "RA---TAN",
      "CTYPE2": "DEC--TAN",
      "CRVAL1": 331.0984070266234,
      "CRVAL2": -56.801257162911156,
      "CRPIX1": 25.0,
      "CRPIX2": 25.0,
      "CDELT1": -0.0017083333333333334,
      "CDELT2": 0.0017083333333333334
    }
  },
  "frames": [
    {
      "index": 0,
      "label": "Эпоха 1",
      "mjd": 60796.1287,
      "date": "2025-05-01",
      "mjd_span": [
        60794.0655,
        60802.2204
      ],
      "det": 3,
      "wave_um": 1.8414,
      "n_exp": 6,
      "image_url": "/cache/27f2642081a832d7/f0.png",
      "flags_url": "/cache/27f2642081a832d7/f0_flags.png",
      "predicted_px": [
        24.996,
        23.383
      ],
      "measured_px": null,
      "measured_candidates": [
        {
          "kind": "largest",
          "xy": [
            14.045,
            8.522
          ],
          "area_px": 26
        },
        {
          "kind": "central",
          "xy": [
            24.915,
            23.421
          ],
          "area_px": 26
        }
      ],
      "obs_ids": [
        "2025W18_1B_0081_2",
        "2025W18_1B_0091_2",
        "2025W18_1B_0175_2",
        "2025W18_1B_0238_4",
        "2025W18_2B_0011_4",
        "2025W19_1B_0344_4"
      ]
    }
  ],
  "trail_url": null,
  "motion": {
    "pm_reference_arcsec_yr": [
      3.982,
      -2.4668
    ],
    "source": "Gaia DR3 / gaiadr3.gaia_source",
    "pm_measured_arcsec_yr": null
  },
  "warnings": [
    "λ в точке цели меняется между эпохами (1.76-2.06 мкм): у объектов с крутым спектром остатки в разностном кадре могут быть не движением",
    "В центре поля есть насыщенные или загрязнённые пиксели (слой флагов)"
  ],
  "credits": "NASA SPHEREx via NASA/IPAC IRSA"
}
```

### 200: job error or empty result

- Job error: `200 {"status": "error", "message": "..."}`. The job is then deleted and a repeated request starts it again.
- No images: `200` with `"frames": []` and a warning that SPHEREx has not yet observed the region. This result is not cached so new SPHEREx releases can be picked up.
- `diff` requires at least 2 epochs; otherwise `frames` is empty and a warning is returned.

## GET /api/export.gif

Animated GIF from completed frames. Parameters are the same as `/api/view`. It uses the cache only: if no result exists, it returns `409 not_ready` (open `/api/view` first).

Response: `image/gif`. Error example: `mock/error_not_ready.json`.

## GET /api/sources

Lists loaded local epoch tables (`tables/` or `SPHEREX_TABLES`). If a table covers the position, epochs are taken from it quickly; otherwise IRSA TAP is queried. A table is an IRSA “Spectral Images” export in `.tbl` format; changed files are re-read without restart.

```json
[
  {
    "file": "spectral-images.tbl",
    "rows": 372,
    "ra": [
      292.32453152743597,
      300.34791431553333
    ],
    "dec": [
      -18.687664822415275,
      -10.772619538696
    ]
  }
]
```

## GET /api/spectrum

On-demand spectrum: photometry at the target position across all exposures from all 6 detectors (`src/processing.build_spectrum`). The same mechanism as `/api/view`: `202` with progress (`stage`: `query`, `download`, `process`, `render`), polling the same URL until `200`, and caching the result.

| Parameter | Description |
|---|---|
| `target` or `ra`+`dec` | As in “Common parameters” |

The target position at each exposure date uses proper motion (SIMBAD); `ra`+`dec` without a name is treated as stationary. Two images per epoch are used for each detector, so the first cold request is long (dozens of cutouts).

Point fields: `wave_um` and `bw_um` (wavelength and bandwidth at target), `flux_ujy` and `err_ujy` (aperture photometry in µJy; `null` if all aperture pixels are bad), `mjd`, `det`, `obs_id`, `ok` (`false` for bad aperture pixels, such as saturation of bright stars; such points are shown but flagged), `flags` (aperture flag names). `n_saturated` is the number of points with `ok = false`. The example includes the first 3 of `n_points`. In addition to `points`, flat arrays `wavelength`, `flux`, `flux_err` contain only reliable (`ok = true`) points for direct frontend plotting; they are omitted from the example.

```json
{
  "status": "ready",
  "points": [
    {
      "wave_um": 0.77603,
      "bw_um": 0.01924,
      "flux_ujy": 166.639,
      "err_ujy": 68.115,
      "mjd": 60968.08653,
      "det": 1,
      "obs_id": "2025W43_1A_0019_1",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    },
    {
      "wave_um": 0.77604,
      "bw_um": 0.01923,
      "flux_ujy": -89.501,
      "err_ujy": 71.447,
      "mjd": 60990.15914,
      "det": 1,
      "obs_id": "2025W46_1A_0287_1",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    },
    {
      "wave_um": 0.85238,
      "bw_um": 0.02151,
      "flux_ujy": 1111.941,
      "err_ujy": 79.613,
      "mjd": 60991.85824,
      "det": 1,
      "obs_id": "2025W46_1A_0811_1",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    }
  ],
  "n_points": 64,
  "n_saturated": 10,
  "credits": "NASA SPHEREx via NASA/IPAC IRSA",
  "bands": {
    "1": [
      0.75,
      1.09
    ],
    "2": [
      1.1,
      1.62
    ],
    "3": [
      1.63,
      2.41
    ],
    "4": [
      2.42,
      3.82
    ],
    "5": [
      3.83,
      4.41
    ],
    "6": [
      4.42,
      5.0
    ]
  },
  "query": {
    "ra": 331.0984070266234,
    "dec": -56.801257162911156,
    "target": "Epsilon Indi B"
  }
}
```

## GET /api/lightcurve

Flux over time in a narrow wavelength window from the same cached spectrum. If the spectrum is still processing, its `202` is returned.

| Parameter | Default | Description |
|---|---|---|
| `target` or `ra`+`dec` | none | As in “Common parameters” |
| `wave_min`, `wave_max` | `±5%` around the centre of detector `det`'s band | window in µm; requires `0 < wave_min < wave_max` |
| `det` | `3` | detector used for the default window |

The window is intentionally narrow: the linear filter produces different fluxes at different wavelengths, and without a window variability would be confused with the spectral shape. As with the spectrum, the response contains flat arrays of reliable points: `dates` (`YYYY-MM-DD`), `flux`, `flux_err`. `warnings` repeats spectrum warnings and appends a warning if no reliable points fall within the window.

```json
{
  "status": "ready",
  "wave_window_um": [
    1.9189999999999998,
    2.121
  ],
  "points": [
    {
      "wave_um": 1.93446,
      "bw_um": 0.04779,
      "flux_ujy": 12360.832,
      "err_ujy": 121.939,
      "mjd": 60836.19048,
      "det": 3,
      "obs_id": "2025W24_1A_0245_4",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    },
    {
      "wave_um": 1.97884,
      "bw_um": 0.04887,
      "flux_ujy": 14348.225,
      "err_ujy": 128.539,
      "mjd": 61164.11762,
      "det": 3,
      "obs_id": "2026W19_1A_0041_1",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    },
    {
      "wave_um": 1.93385,
      "bw_um": 0.04716,
      "flux_ujy": 12295.82,
      "err_ujy": 118.852,
      "mjd": 61311.33729,
      "det": 3,
      "obs_id": "2026W40_1A_0070_2",
      "ok": true,
      "flags": [
        "SOURCE"
      ]
    }
  ],
  "query": {
    "ra": 331.0984070266234,
    "dec": -56.801257162911156,
    "target": "Epsilon Indi B"
  }
}
```

## Static files

- `/cache/<key>/f<i>.png` and `/cache/<key>/f<i>_flags.png`: images from `image_url` and `flags_url`.
- `/`: frontend. It accepts URL parameters, e.g. `/?target=barnard&fov=3&det=3&mode=coadd`.

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `CACHE_DIR` | `cache/` | Cache directory (results, FITS cutouts, reference data) |
| `PREWARM` | `1` | `0` disables warming the five presets at startup |
| `CENTER_MJD` | today | Fixes the centre recalculation date (stable cache keys for demos) |
| `SPHEREX_TABLES` | `tables/` | Directory of local epoch tables |

## Demo cache

To avoid live downloads during a demonstration, the cache is prepared in advance with `scripts/prepare_demo_cache.py` (all presets, `coadd` and `diff`, spectrum and light curve), and the server starts with `start_demo.sh`:

```
CENTER_MJD=61358 PORT=8000 bash start_demo.sh <папка_кэша> 61358
.venv/bin/python scripts/prepare_demo_cache.py --base http://127.0.0.1:8000
```

The cache key depends on field centre, which is recalculated for `CENTER_MJD` (default “today”). The prepared cache works only with the same date: the ready cache uses `CENTER_MJD=61358` (November 14, 2026). Without a fixed date the keys change daily.

## Known limitations

- Proper motion from measured centres (`motion.pm_measured_arcsec_yr`) is not calculated.
- Spectrum and light curve use the `src/Exposure` loader (its own concurrency limit and 503 retries), not `api/irsa.py` search; local epoch tables do not affect them.
- The methane index is not implemented.
- Very bright targets saturate images: all Barnard spectrum points have saturation flags (flux is measured for only 8 of 48 points); R Dra has no flux in any point. These points are returned but marked `ok: false`.
- The default light-curve window is narrow (±5% around the detector midpoint), so some targets have no points in it (cached: RW Aur, R Dra and Barnard have 0, WISE 0855 and eps Ind B have 3 each). Set `wave_min` and `wave_max` manually to see more.
- A cold request takes tens of seconds to several minutes (about 4 MB and 13 seconds per cutout); cached responses take milliseconds.
- The object centre uses a 3σ threshold and first moment, without PSF fitting or neural networks; the threshold may need adjustment on very noisy frames (constants in `api/centroid.py`). `err_px` includes noise only.
- Touching objects are not separated (`elongated`); an object near a mask edge has a shifted centre (`touches_nodata`).
- In `diff`, a nearby negative lobe eats into the positive lobe: at ~5 pixels displacement the centroid bias is about 0.35 pixels, and greater for smaller displacements.
- Only constant residual background is subtracted from each exposure; gradients remain in large fields.

## Dilayla's contract clarification (November 14, 2026)

`/api/resolve` retains `ra`, `dec`, PM and `source`, and adds `source_id` (string), `source_url`, `ref_epoch_format`, `ref_epoch_scale`, `coordinate_model`, `assumptions`, and, for the team catalogue, component, limitations and available errors. `ref_epoch` retains the original value: for WISE this is **57633.39 MJD**, not a year. Unknown values return `null`; conditional RV=0 remains only a model assumption. WISE's unspecified source epoch scale remains `null`; UTC use appears in `assumptions`.

The five targets use `sky_motion.predict_position`: geometric geocentric direction in ICRS axes, without aberration or spacecraft-orbit correction. This is a prediction, not a measurement. For arbitrary Sesame responses, J2000 is not treated as a verified motion reference epoch: coordinates are fixed (`fixed_catalogue`), and PM is retained as a catalogue value. With a reliable epoch and PM the adapter supports motion propagation using the specified scale. The full name is encoded in the outgoing URL; service errors are distinguished from a missing name.

`/api/view` adds `target_info`, `center_mjd` (MJD UTC for the centre of a name query), `data_status`, and `exposure_counts` (`selected`, `failed`, `used`; absent before exposure search). States: `no_archive_results`, `download_failed`, `unusable_data`, `partial`, `single_epoch`, `available`. For diff with one epoch there are no frames and a warning explains the limitation. An empty archive response does not prove absence of observations for all time. A poor-quality measurement does not populate `measured_px`; candidates are retained separately. An out-of-frame prediction is hidden with a warning. Pixels are zero-based `[x,y]`, y downward.

FITS time is normalized to MJD UTC accounting for `TIMESYS`; `MJD-AVG` takes priority, followed by `DATE-AVG`. A frame contains exposure `time_sources` with source keys and assumptions. `archive_hint` may be an observation's lower bound. A missing TIMESYS means an explicit UTC assumption. Unknown time must not be interpreted as a verified midpoint.

Target information and spectrum warnings are retained when obtaining a light curve. Cache version 11 invalidates previous results; fixed and tracked queries with the same centre have different keys. Original FITS are retained. On the main page, plots start via buttons; an empty form loads only presets. URL parameters include `max_per_epoch`, `api`, `mock`. Synthetic data are explicitly labelled.

Reproducible limited check:

```sh
PREWARM=0 .venv/bin/python -m uvicorn api.main:app --host 127.0.0.1 --port 8011
.venv/bin/python scripts/validate_targets.py --images --timeout 120
```

Validation responses are in `data/validation/responses/`; a real HTTP request may use the server cache. `mock/` examples use a historical format version.

### Supplement: coverage, quality and freshness

Runtime cache version is **11**. Each frame contains:

- `quality`: `status` (`limited`/`unchecked`), `verified=false`, finite-pixel fraction `coverage_fraction`, and `reasons` codes;
- `prediction_status`: `in_frame`, `outside_field`, or `no_model`;
- `wavelength_samples`: actual λ/bandwidth for every used exposure, `obs_id`, and source `spectral_wcs` or `approximate_detector_position`;
- `time_sources`: only exposures actually used in the frame, UTC normalization, and time assumptions;
- overall `exposure_counts`: `selected`, `failed`, `used`, `rejected`;
- overall `spectral_sampling`: `wave_min_um`, `wave_max_um`, `channel_width_um` (median channel width), `n_samples` (number of usable λ samples of used exposures), and `different_channels`, `true` if target wavelengths in different frames differ by more than a channel width (with a `warnings` notice that brightness differences may depend on spectral shape and variability is not confirmed).

For `/api/view`, wavelengths and channel width in `spectral_sampling` are rounded to 4 decimal places. Without usable samples `n_samples=0` and the other four fields are `null` (including `different_channels`). When channel width is unknown, comparison uses 5% of mean λ.

For `/api/lightcurve`, diagnostics concern only usable points in the selected window: `wave_min_um`, `wave_max_um`, and `channel_width_um` are unrounded; `n_samples` is absent. `different_channels` is boolean; with no usable points it is `false` and the three numeric fields are `null`. If channel width is unknown, `different_channels=false` does not establish matching bands. For either endpoint, `false` does not establish suitability for quantitative flux comparison.

With `measurement.target_association=prediction_outside_field`, the model places the target outside the frame: `measured_px=null`, `offset_from_prediction_px=null`; visible sources remain in `measured_candidates`, and the interface explicitly reports no prediction in the field.

`quality.verified=false` remains even with complete coverage: this is not scientific certification of an image. `measured_candidates` and an unambiguously algorithm-selected centroid do not establish the star's identity. The difference scale is `robust_residual_scale`, without calibrated probability/significance. A difference with no common usable pixels is not presented as a useful PNG frame.

`processing.registration` describes only header-WCS reprojection: `name=header_wcs_reprojection`, `independently_verified=false`, `astrometric_precision_arcsec=null`. Independent registration using reference stars and measurement accuracy are not confirmed. Warnings reach the main page: centroids are candidates, and brightness changes and difference residuals require comparable bands/responses and quality checks. This applies even when λ differs by less than 5% or mean λ matches. WCS-WAVE provides a band for viewing; quantitative analysis requires checking spectral calibration/response and PSF.

Spectral points add `quality_reasons`, `coverage`, `wavelength_source`, `time_source`, and `time_assumptions`. Incomplete apertures/background annuli and aperture contamination produce `ok=false`. The cutout includes a margin for motion and the annulus on selected dates. Unavailable spectral calibration gives `wave_um=null`, not JSON NaN. `n_saturated` counts OVERFLOW/NONLINEAR flags; `n_flagged` counts all unusable points. Flux error includes aperture-pixel variance and background-annulus median error (`method.error_components`); aperture correction, flux calibration, target-position error, and channel shape are excluded (`method.error_not_included`). A background annulus more than half masked yields `background_annulus_mostly_masked`. The light curve returns `spectral_sampling` and warns if points use different channels. `different_channels=false` is only a centre/width heuristic result, not confirmation of identical responses or usable photometry.

Normal search uses the archive. `SPHEREX_SOURCE=tables` explicitly enables a limited local set with keys separate from archive keys. Completed results and exposure lists refresh after 3600 seconds, partial and empty results after 60 seconds; settings are in [DEPLOYMENT.md](DEPLOYMENT.md).

HTTP 503 `busy` means the queue is full; `timeout` means the job time limit was exceeded. The same job is not restarted while already-started network operations finish. The client may retry after an error. Two identical requests use one job; rechecking the cache under a lock eliminates the result-publication race.

### TAP failure and SIA2 fallback

`/api/view` first sends one GET request to IRSA TAP with a 12-second limit. On 502, timeout or unusable CSV, official IRSA SIA2 search automatically queries collections `spherex_dr1` and `spherex_qr3`. Duplicate observations on the same detector are combined with DR1 priority. This is a live archive request, not an old result. Direct SIA2 FITS links are downloaded without incorrectly attempting to parse them as XML datalink.

`archive_source` in the view response is `tap`, `sia2`, or `tables`. With `sia2`, `warnings` notes that new observations may appear in SIA2 roughly a day later than TAP. Source and warning are cached even for an empty set. `QUERY_STATUS=ERROR` and truncated `OVERFLOW` responses are not represented as absent coverage or cached as successful results.

One failed TAP probe disables TAP for all jobs in the process for 300 seconds; TAP is then checked again. `TAP_PROBE_TIMEOUT` and `TAP_RETRY_SECONDS` configure these intervals. FITS downloads retain retries for 408, 429, and 5xx; permanent 4xx responses are not retried.

Spectrum, light curve and Python library use SIA2 directly. Their searches also have a finite timeout (45 seconds) and up to three attempts for temporary errors, with 2- and 4-second delays. Search failure is classified as archive failure, not a local processing error. FITS downloads retain their own retries and cache.

References: [IRSA TAP documentation](https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html) and [official SPHEREx data access](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/).

## Rate limiting and public errors

`/api/view`, `/api/spectrum`, `/api/lightcurve`, `/api/resolve`, and `/api/export.gif` share a token bucket per network client. By default, bursts of up to 120 requests are allowed, then the budget refills at 120 per minute. `API_RATE_LIMIT_PER_MINUTE` sets a positive limit; 202 polling counts toward the same budget. PNGs, the home page, healthcheck and CORS OPTIONS are exempt.

When exhausted: HTTP 429, `error.code=rate_limited`, and `Retry-After` in seconds; processing does not start. 503 `busy`/`timeout` responses also include `Retry-After: 5`. CORS allows clients to read this header.

Unexpected HTTP 500 and background-worker errors contain a generic message. Exception details remain in server logs and are not sent to visitors. The limiter uses `request.client` and does not trust arbitrary X-Forwarded-For. Behind a reverse proxy, trusted forwarding addresses must be configured explicitly on the server; otherwise clients of one proxy share a budget. Limiter state is stored in one process and capped at 4096 client entries; multiple replicas require a shared ingress limiter.

SOURCE alone is not removed even by the isolated-outlier heuristic: without independent verification, such a pixel is not declared defective. An accompanying HARD_BAD flag still excludes the pixel.
