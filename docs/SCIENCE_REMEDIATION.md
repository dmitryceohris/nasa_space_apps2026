# SPHEREx remediation: G1–G7

Date: November 15, 2026. Changes were made locally; no publication or GitHub submission was performed.

## Goals and acceptance

| Goal | Implemented | Verification |
| --- | --- | --- |
| G1 — masks and selection | PHANMISS, DICHROIC, SATELLITE_HALO; SOURCE retained at the target and excluded from the background; incomplete apertures are not reconstructed; ok=False is excluded from science arrays | Regression tests for masks, contamination, incomplete regions, and missing FLAGS |
| G2 — observations and calibrations | Unified registry with source sample membership, release, version, URL, SHA-256 and cutout files; CWAVE/CBAND and area maps for all six detectors | Cutout coordinates checked against the full detector; real FITS |
| G3 — validation | Shared validated/provisional/rejected statuses and reasons; separate photometry and astrometry checks; diagnostic points behind an explicit switch | Python, API, interface and FITS export |
| G4 — ePSF and spectral selection | Native-pixel fit, local background, VARIANCE, masks, pixel area, joint fitting of neighbours, errors and residuals | Independent artificial stars, different subpixel positions, blending, 40 noise realisations |
| G5 — independent astrometry | Full detection list without selection by prediction, free ePSF centroids, registration against reference stars and covariance propagation | End-to-end test of the target and eight independent reference stars; outliers and insufficient references |
| G6 — model and statistics | Observer geometry, covariance between epochs, fixed-parallax error, free-fit geometry checks, correlated flux errors | Recovery of prescribed motion, degenerate geometry and invalid covariances |
| G7 — real data | Five primary targets, all six detectors, original Barnard samples and three additional stars | Reproducible report with accepted, rejected and unavailable inputs |

Implementation and algorithm checks are complete. The scientific fitness of each real result is determined by its own status: code readiness does not confirm a measurement of motion or variability.

## Results on real data

The [measurement report](../data/validation/science-remediation/measurements.json) contains **245 inputs/attempts**, including repeated cutouts of the same products: **225 processed FITS**, **20 unavailable**, **10 accepted WISE 0855 photometric results**, and **215 rejected results**. Repeated cutouts are not independent epochs. The registry tracks physical products separately from files and sample membership.

For the original three D1/DR1 Barnard observations, calibrated values were obtained at the geocentrically predicted position:

| obs_id | λ, µm | Width, µm | Result |
| --- | ---: | ---: | --- |
| 2025W35_1A_0432_4 | 0.749428 | 0.017463 | Damaged core / poor fit |
| 2025W40_1A_0471_1 | 1.031579 | 0.027205 | Damaged core / poor fit |
| 2026W11_1A_0539_2 | 0.925587 | 0.023908 | Damaged core / poor fit |

These bands differ substantially; the three points cannot be combined into one light curve. The original sample files and physics_observations.csv are preserved; their contents are explicitly separated in the [registry](../data/validation/science-remediation/observations.json).

The checked sample contains **no confirmed light curve or proper motion** for any target. Masks and saturation matter for Barnard, R Dra and RW Aur; model quality and the unresolved system matter for Epsilon Indi B. WISE has usable forced photometry, but this is not an independent detection of motion. This conclusion concerns the checked sample and does not prove that no suitable frames exist elsewhere in the archive.

The search covers metadata from all six detectors; pixel-level checks are limited to six distinct dates per detector and additional original samples. Search responses combine newly saved responses from the current run and saved November 14 responses used when the archive failed. SHA-256 hashes of raw responses are retained; current archive completeness is not verified.

Five ESA Gaia reference-star queries timed out. They are recorded as unavailable; registration of real frames is not claimed to be confirmed. Saved frames for three additional stars were checked, without attributing a new independent astrometric solution to them.

## Rules and limitations

- Preliminary aperture photometry is always provisional; incomplete or contaminated apertures are rejected without scaling the remaining pixels.
- The main method automatically uses ePSF when calibrations are available. Fitting accounts for pixel response once: the ePSF is already a pixel-convolved model. Normalisation is determined from the full model before masking or cropping.
- Flux uncertainty contains a statistical component and 3% of the flux in quadrature. **3% is an explicitly documented project assumption**, not a universal IRSA standard established by this work. Calibration correlation between epochs is not guessed.
- Detection requires S/N≥5 using the statistical error. Model checks require full rank, no active bounds, χ² p≥0.001, and a mean absolute normalised residual within 2.5 px of ≤2. These are project thresholds, tested on synthetic data but not tuned to real targets.
- Light curves require full spectral responses with provenance. CWAVE/CBAND overlap <90% rejects a comparison; proximity of intervals does not validate it. For supplied full responses, integrated overlap ≥98% is required pairwise within each group. These are project tolerances. Without responses only diagnostic points are available.
- `constant_flux` requires accepted comparable photometry and an explicit covariance or calibration correlation model. Duplicate physical observations are rejected.
- `GEOCENTER` defines the spacecraft coordinate origin but does not establish axis orientation. Without documented orientation, the geocentric model remains provisional and the maximum orbital correction is estimated. In the current sample its upper bound reaches approximately 0.000026″.
- Free parallax remains exploratory without six distinct dates, a one-year span, three annual phase quarters, or matrix stability. Catalogue parallax includes propagation of its published error.
- `Motion.validation` separates an exploratory fit from accepted motion. Old common-grid centroids and assisted measurements do not become independent astrometry. `Motion.pm`/`to_numpy()` values remain diagnostic when provisional.
- For Epsilon Indi B and RW Aur, a single unresolved model does not validate component astrometry. Flux in the report refers to the system.
- Image differences and FITS export use `robust_residual_scale`, not statistical σ. The VARIANCE extension explicitly refers to the coadd, not the normalised difference.

## Reproduction

```sh
# Network: prepare pinned calibrations and, if ESA is available, reference stars.
uv run --locked python scripts/prepare_calibrations.py
uv run --locked python scripts/prepare_references.py

# Offline: repeat measurements using the local FITS cache and saved searches.
uv run --offline --locked python scripts/accept_science.py --per-detector 6

# Separately authorised refresh of the search and missing cutouts.
uv run --locked python scripts/accept_science.py --refresh-archive --per-detector 6 --download

# Regression tests and end-to-end API replay of saved FITS.
uv run --offline --locked pytest -q
uv run --offline --locked python scripts/accept_processing.py \
  --fits-cache cache --output data/validation/science-remediation/api-replay.json
git diff --check
```

On another computer, restore FITS files using the URLs and SHA-256 hashes in the manifests. Calibration files are in `cache/calibrations`; their provenance and hashes are also recorded in the [manifest](../data/validation/science-remediation/calibrations.json). Set `SPHEREX_CALIBRATION_MANIFEST` for another location. Calibration matching is explicitly pinned by release and processing version; it is not inferred from observation date.

For local FITS, pass `product_url=...` when loading to preserve the release and version. `spx.native_track(exps, target, reference_catalog=provider)` accepts a catalogue propagated to each frame time and its covariances; prepared catalogues from `cache/references` are loaded automatically. `spx.spectrum(...).to_numpy()` contains only accepted photometry; `only_ok=False` enables diagnostics. For light curves, request diagnostics with `include_diagnostics=True`.

The `/api/spectrum` and `/api/lightcurve` routes are retained. Added fields include `validation`, `quality_counts`, `diagnostic_points`, `series`, provenance and error components. The result cache version was increased to 12. Science arrays in the old cache are recalculated.

Final test result: **314 passed, 1 skipped**. The skipped test is an existing live archive check. End-to-end API replay: **47 real inputs, 36 scenarios, 188 frames**, without archive requests. The [replay result](../data/validation/science-remediation/api-replay.json) contains hashes of the tested code. Browser verification uses controlled synthetic responses and the real interface: an accepted point, explicit enabling of diagnostics, prevention of connecting incomparable bands, and a width of 375 px. It does not confirm availability of external archives.

Implementation references: [IRSA products and calibrations](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/), [ePSF and cutout coordinates](https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/), [photometry and FLAGS](https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/sp.html).
