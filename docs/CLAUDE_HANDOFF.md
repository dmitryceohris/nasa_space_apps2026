# Claude handoff: five-target SPHEREx physics

Branch: `codex/physics-five-targets`

## Delivered

- Added `sky_motion.py`, a shared five-target predictor that reuses the existing
  `coordinate_at` and time handling from `barnard_motion.py`. It adds geometric
  geocentric annual parallax with Astropy's built-in Earth ephemeris.
- Prepared predictions for 15 real SPHEREx observations: three epochs for each
  challenge target. `output/predictions.csv` is reproducible offline from
  `data/physics_observations.csv`.
- Saved SIA search results, selected observation metadata, cutout URLs, FITS
  headers, local FLAGS evidence, retry outcomes, and selection rationale under
  `data/physics_evidence/`.
- Added `docs/PHYSICS.md` with model assumptions, all five scientific
  expectations, validation status, the proposed CH4-K diagnostic, and slide text.
- Updated the project README.

## Verification

`uv run --offline --locked pytest -q` passed: 57 tests, including all 38 existing
Barnard tests and tests for the generalized calculation, parallax signs, epoch
selection, footprint checks, FITS flags, and offline CSV generation.

The original Barnard dates reproduce 5.594198 arcsec barycentrically and
5.454815 arcsec with geocentric geometric parallax. The five-target CSV has
three rows per target, with every selected position more than 10 arcsec inside
the searched footprint.

## Archive follow-up

Of the 15 chosen frames, 3 passed the local fast FLAGS check, 4 were flagged,
and 8 remain unchecked after IRSA returned HTTP 503 on both bounded attempts.
Their positions are predictions, not measurements. Do not present the proposed
methane index as a measured value; the required photometry has not been done.
See `docs/PHYSICS.md` and `data/physics_evidence/selection.json` for details.

## Reproduce

```sh
uv run --offline --locked python sky_motion.py --observations data/physics_observations.csv --output output/predictions.csv
uv run --offline --locked pytest -q
```
