# GIFs and statistics for five SPHEREx targets

Saved real FITS observations reprocessed using the application code on 15 November 2026.
No new IRSA archive search was performed. All GIFs were generated using `/api/view` and `/api/export.gif`.

## Application view

| Target | GIF frames | Exposures | Detector | Original exposure dates, UTC | Baseline, days |
|---|---:|---:|---|---|---:|
| Barnard's Star | 5 | 23 | D3 | 2025-08-15 — 2026-08-20 | 370.0 |
| Epsilon Indi B | 7 | 32 | D3 | 2025-04-29 — 2026-09-28 | 517.3 |
| WISE 0855 | 4 | 18 | D6 | 2025-05-05 — 2026-05-11 | 370.4 |
| RW Aurigae | 3 | 18 | D2 | 2025-09-08 — 2026-09-15 | 372.1 |
| R Draconis | 5 | 30 | D2 | 2025-06-08 — 2026-06-23 | 379.7 |
| **Total** | **24** | **121** | | | |

The frame count is the number of epochs after grouping. A frame can combine multiple exposures.
The date in a GIF is the mean frame date; the table lists the dates of the earliest and latest original exposures.
The ring shows the catalogue model prediction. Each GIF uses a fixed field across epochs.
A consistent brightness scale is used throughout each sequence, with no interpolation between dates.
Violet pixels indicate no usable data. Magnification does not add resolution.
Different spectral channels and masks can affect appearance. Validated astrometric frames: 0.

## Separate scientific assessment

This assessment includes different exposures and all detectors; it does not cover the same sample as the GIFs.
Source: `data/validation/science-remediation/measurements.json`; checked at `2026-11-15T06:32:14.701` UTC.
Repeated cutouts of the same product were excluded using the full product URL.

| Target | Unique products | Accepted photometry | Rejected | Unavailable |
|---|---:|---:|---:|---:|
| Barnard's Star | 42 | 0 | 39 | 3 |
| Epsilon Indi B | 45 | 0 | 40 | 5 |
| WISE 0855 | 41 | 6 | 28 | 7 |
| RW Aurigae | 38 | 0 | 37 | 1 |
| R Draconis | 48 | 0 | 44 | 4 |
| **Total** | **214** | **6** | **188** | **20** |

Accepted WISE 0855 photometry was measured at the catalogue position.
Neither proper motion nor variability is confirmed for any of the five targets in this dataset.
These numbers describe limited saved samples, not the entire SPHEREx archive.

## Exact WISE 0855 coordinate input check

RA `08h55m10.831s` and Dec `-7d14m42.53s` were passed to `/api/view` without a target name.
The additional `wise0855-coordinates/blink.gif` sequence contains 4 frames,
4 exposures, D6, and a 6′ field. This is a fixed field without a motion model.
It is excluded from the total of 24 frames and 121 exposures in the five main GIFs.

## Files

- `index.html`: the gallery, available at `/challenge-demo/` when the application is running.
- `<target>/blink.gif`: version with title, date, wavelength and catalogue prediction.
- `<target>/api.gif`: original application GIF export without labels or enlargement.
- `<target>/api_response.json`: complete application response; its `/cache/` URLs refer to the isolated run.
- `<target>/frame_*.png`: local copies of frames and masks for viewing without the temporary cache.
- `<target>/inputs.json`: source products and checksums of the FITS files and saved metadata.
- `statistics.json`, `statistics.csv`: application statistics and the separate scientific assessment.

To reproduce from the project root, run `.venv/bin/python scripts/export_challenge_gifs.py`.
The original local FITS cache listed in the manifests is required. The export itself disables network access.
