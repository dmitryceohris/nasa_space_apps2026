# Physics of the five SPHEREx targets

Prepared November 14, 2026. These are predictions from catalogue astrometry, not measured centroids or fluxes, nor evidence that the objects were detected.

## Running from the repository root

Retrieve metadata and small cutouts **only for quality checks**:

```sh
uv run --locked python sky_motion.py --prepare-observations --output data/physics_observations.csv
```

Recalculate offline using saved times:

```sh
uv run --offline --locked python sky_motion.py --observations data/physics_observations.csv --output output/predictions.csv
uv run --offline --locked pytest -q
```

No new dependencies are required. Input parameters are in `OBJECT_PARAMETERS.csv`; their provenance and component identities are discussed in `OBJECT_PARAMETERS.md`. The existing `barnard_motion.py`, its interface, notebook and results are retained.

To call from other code:

```python
from astropy.time import Time
from sky_motion import predict_position

p = predict_position("wise_0855", Time("2026-04-30T00:00:00", scale="utc"))
ra, dec = p.geocentric.ra.deg, p.geocentric.dec.deg
print(p.assumptions)
```

Accepted short names: `barnard`, `epsilon_indi_b`, `wise_0855`, `rw_aur`, `r_dra`; exact names from the catalogue's `object` column are also accepted.

## Model and important limitations

`coordinate_at` from `barnard_motion.py` propagates proper motion using `SkyCoord.apply_space_motion`, including the perspective effect when RV is supplied. Annual parallax is calculated separately as **r_geo = r_star,bary − r_Earth,bary**; Earth's position comes from `get_body_barycentric("earth", time, ephemeris="builtin")`. RA/Dec give this vector's direction in axes parallel to ICRS. [Astropy: space motion](https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html), [ephemerides](https://docs.astropy.org/en/stable/coordinates/solarsystem.html).

This is a **geometric geocentric direction**, not a full apparent GCRS position: aberration, gravitational deflection, the SPHEREx orbit, and additional light-travel-time corrections are not added. Check the specific WCS conventions before precise comparison with measurements; aberration already corrected in the registration must not be added again.

- For four objects without RV, the conditional model `RV=0` is used with the `rv_assumed_zero` flag. This is neither a velocity measurement nor an estimate of its uncertainty.
- WISE's reference epoch MJD 57633.39 is provisionally treated as UTC: `reference_epoch_scale_assumed_utc`. Unknown correlations and RA error convention are not turned into a fictitiously precise covariance.
- Epsilon Indi B coordinates describe the Gaia Ba source, not the Ba/Bb barycentre or infrared photocentre. RW Aur uses the eastern component called A in the selected literature. Elevated RUWE is noted in `assumptions`.
- The new calculation does not produce formal errors; extending covariance propagation is deferred. The existing Barnard error calculation continues to work separately.

Orders of magnitude of omitted corrections, without a separate orbital implementation:

- At a SPHEREx distance of about 7000 km from the geocentre, the correction for Barnard is at most on the order of `0.547″ × 7000 / 149597870.7 = 0.000026″`, or about **0.00003″**, equivalent to **0.000004 pixels** at 6.15″/pixel.
- For WISE, a time difference of roughly one minute gives `8.152″/year × 60 / 31557600 ≈ 0.000016″`, about **0.00002″**.

Compared with a *hypothetical*, not yet measured centroid error of 0.5″, both effects are about four orders of magnitude smaller. Earth's geometric **annual** parallax reaches tenths of an arcsecond and is therefore included.

## Epochs, quality and reproducibility

SIA DR1/QR3 responses are saved in `data/physics_evidence/`; for Barnard the already saved response in `data/irsa_barnard_20261114/` was used. Absence of `QUERY_STATUS=OVERFLOW` was checked. A 0.05° search radius around the catalogue position covers the targets' motion over the available interval.

Epochs are separated by gaps greater than 30 days. One frame is selected per cluster, with at most three clusters: first, middle and last. For one continuous cluster, three equal time windows are used; these are not separate surveys. Candidates prefer DR1, then D6 for WISE and D1 for the others; saved Barnard frames take priority within the corresponding epoch.

Coverage is checked using spherical `s_region` at the propagated position, with a required margin greater than 10″. All 15 selected frames have margins greater than 100″. Time is taken from `MJD-AVG`/`TIMESYS`, then `DATE-AVG`. For SIA the interval midpoint is calculated on the continuous TAI scale; the result is stored as MJD UTC and receives `sia_time_scale_assumed_utc` because the scale has not been independently confirmed. Saved headers contain `MJD-AVG` and `TIMESYS=UTC`.

A 5×5 FLAGS region near the predicted position is checked for `OVERFLOW`, `NONLINEAR`, `SUR_ERROR`, `NONFUNC`, `MISSING_DATA`, `HOT`, and `COLD`; bits are read from the product's `MSKNnnnn`. IMAGE values are also checked for finiteness. On failure one additional candidate is allowed, preferably from another detector. [IRSA: FLAGS and spectral calibration](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/), [small-cutout service](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/).

Status meanings:

- `usable`: this quick check passed; detection, photometric suitability, and absence of other contamination are not established.
- `flagged`: bad pixels were found in the checked region. This does not necessarily establish saturation of the star itself; contamination may be nearby.
- `unchecked`: the check was not completed, for example after two HTTP 503 responses.
- `no_coverage`: the predicted position is outside the allowed coverage region.

`selection.json` records selected frames and both attempts, URLs, service errors, local FLAGS values, bit descriptions, and pixel coordinates. IMAGE/FLAGS headers for successful downloads are saved alongside it. Full images were not downloaded; the small cutouts obtained are not retained. They may contain a shared ePSF table, so a response of about 4 MB is possible even with a small IMAGE array. Timeout is 15 seconds, with at most one retry per request. The total network stage is limited to 15 minutes; a partial result is saved after each target.

### Saved result

There are **15 observations: 3 usable, 4 flagged, 8 unchecked**. The eight unchecked frames were retained because of IRSA HTTP 503 responses, with explicit SIA time assumptions.

| Target | Selected epoch dates UTC | Checks in date order | Predicted geocentric displacement, first → last |
| --- | --- | --- | ---: |
| Barnard | 2025-08-27; 2026-03-11; 2026-09-02 | flagged; unchecked; flagged | 10.549699″ ≈ 1.7154 pixels |
| Epsilon Indi B | 2025-05-20; 2026-04-27; 2026-09-27 | usable; usable; unchecked | 6.162287″ ≈ 1.0020 pixels |
| WISE 0855 | 2025-05-05; 2025-11-14; 2026-04-30 | unchecked; unchecked; unchecked | 8.031165″ ≈ 1.3059 pixels |
| RW Aur | 2025-09-04; 2026-03-04; 2026-09-09 | usable; flagged; flagged | 0.025383″ ≈ 0.0041 pixels |
| R Dra | 2025-06-06; 2025-11-20; 2026-05-13 | unchecked; unchecked; unchecked | 0.019594″ ≈ 0.0032 pixels |

For Epsilon Indi B, the frames passing the quick check are from **different detectors D6 and D4**. Their coordinates may be compared after checking registration, but their flux difference does not test variability in a single band.

Check on the **previous three Barnard dates**: 2025-08-27, 2025-10-01, 2026-03-11. The prediction without annual parallax is **5.594198″**, and with geocentric parallax **5.454815″**. The east component changes from **−0.431451″** to **+0.603501″**, and the north component from **+5.577535″** to **+5.421328″**. These numbers refer to the original dates, not the new set of three.

### `output/predictions.csv` columns

`target, epoch_id, obs_id, detector, release` link the row to a frame; `mjd` is MJD UTC, `time_source` is `fits_mjd_avg`, `fits_date_avg`, or `sia_midpoint`. `ra_bary, dec_bary, ra_geo, dec_geo` are degrees. `d_east_arcsec, d_north_arcsec` are spherical offsets relative to **the first selected epoch of the same target**, not the difference between geocentric and barycentric positions. `shift_px` is angular distance from the first epoch / 6.15″, an approximate scale rather than a measured motion in pixels. `coverage_margin_arcsec`, `quality_status`, and `assumptions` describe coverage, checks and assumptions. An empty value means unavailable data, not zero. Frame links are in `data/physics_observations.csv` and `selection.json`.

## What to expect from the five objects

The motion and parallax values below are calculated or read from one source row per object; they are not new SPHEREx measurements.

| Target | μ, ″/year | Catalogue parallax, ″ | Scientific expectation and limitation |
| --- | ---: | ---: | --- |
| Barnard | 10.39335 | 0.54698 | Motion mostly northward. Parallax also changes the direction of inter-epoch displacement. Check bad pixels before measuring the centroid. |
| Epsilon Indi B | 4.68416 | 0.27066 | Southeastward motion; Ba/Bb are T1–1.5 and T6 dwarfs. They are not reliably resolved in these images. Gaia's Ba model does not give the pair's exact IR photocentre. |
| WISE 0855 | 8.15155 | 0.43900 | Motion mostly westward. Signal in a particular channel is not guaranteed; preferring D6 is not a detection. |
| RW Aur | 0.02506 | 0.00545 | Motion over the selected interval is small; look for flux changes in comparable bands. Account for blending and ambiguous A/B designations. |
| R Dra | 0.02154 | 0.00125 | The main effect is Mira variability, not detectable motion. Do not present the literature period of about 246.12 days (AAVSO) as measured from three epochs. |

Astrometry sources: four Gaia DR3 rows and [Kirkpatrick et al. (2021), Table 5](https://arxiv.org/html/2011.11616#S4.T5) for WISE. Epsilon Indi Ba/Bb: [Chen et al. (2022)](https://arxiv.org/html/2205.08077). RW Aur variability: [Koutoulaki et al. (2019)](https://arxiv.org/abs/1903.06423). R Dra classification and literature: [AAVSO](https://targettool.aavso.org/?ac=on&cv=on&ev=on&exo=on&het=on&lpv=on&settype=true&spv=on&yso=on).

## Methane: proposed method, no numerical result yet

The primary target is Epsilon Indi B, but the observed spectrum combines Ba/Bb rather than providing separate indices for each component. The proposed literature **CH₄-K** index lies entirely within D3:

**I = ∫[2.215, 2.255] Fλ dλ / ∫[2.080, 2.120] Fλ dλ**, wavelength in µm.

The numerator is the absorption region, the denominator a comparison region. The index is dimensionless; a lower value indicates deeper relative absorption. Definition: [Burgasser et al. (2006), Table 3](https://arxiv.org/html/astro-ph/0510090#S4.T3).

For SPHEREx this is an **approximate diagnostic**, not a ready classification: at R≈41, channel width around 2.2/41 ≈ 0.054 µm exceeds each 0.040 µm window. The literature integrals cannot automatically be recovered from two channel-centre values. Actual spectral responses and a check that measurements cover the windows are needed. Preferring one detector reduces cross-calibration problems but does not fix insufficient spectral resolution. [SPHEREx bands and resolution](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/).

Before calculation, the team must extract background-subtracted fluxes under comparable conditions and account for flags, errors and responses. MJy/sr in IMAGE is surface brightness Fν, not a ready Fλ for an individual star; photometry, units, solid angle and correct Fν↔Fλ conversion are required. Do not show a numerical index if coverage or signal is insufficient. For an ordinary ratio require a positive denominator with S/N≥5 and numerator with S/N≥3; otherwise report “insufficient data” rather than inserting zero. These are working demonstration thresholds, not a universal methane-detection scale.

Additional CH₄-H windows 1.635–1.675 / 1.560–1.600 µm lie on D3/D2, so the index also depends on cross-calibration and is approximate. For WISE, the region around **3.3 µm within D4** is a possible diagnostic only with sufficient signal. The 4.5/3.3 µm colour is not an unambiguous detection or measurement of methane abundance. Fluxes and index values **were not extracted** in this work.

## Text for the science slide

**The source position is a testable prediction.** We propagated catalogue astrometry for all five targets to actual SPHEREx times and separately added geometric annual parallax. Over the original 196.354-day Barnard interval, expected displacement changes from 5.594″ to 5.455″; its direction changes as well. For the new selection, 15 predictions and local quality checks are saved: 3 frames passed the quick check, 4 have bad flags, and 8 remain unchecked because of archive 503 responses. For RW Aur and R Dra, the main expected effect is a brightness change that must be checked in comparable bands. **Neither a predicted marker nor a literature methane index replaces measurement from the data.**

Measured displacements/fluxes, their errors, and agreement with the prediction are to be added after independent image processing by the team.
