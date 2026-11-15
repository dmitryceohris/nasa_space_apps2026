# Calculating the position of Barnard’s Star

The script `barnard_motion.py` and notebook `output/jupyter-notebook/barnard_motion.ipynb` use the same implementation. The source is the existing **Gaia DR3 4472832130942575872** row in `OBJECT_PARAMETERS.csv`. The parameter table is unchanged. The result is **barycentric astrometric ICRS coordinates**, without annual parallax or a correction for the SPHEREx observer.

Update of November 14, 2026: three IRSA observations were found separately, their FITS headers checked, and positions calculated. An overview of the results is in [README.md](README.md). The “Provided files” section below describes the initial review of the supplied export.

## Installation and execution

Open a terminal in the project directory:

```sh
cd nasa_space_apps2026
uv sync --locked --python 3.12.13
```

A separate `.venv` environment and installed Python **3.12.13** are used. Direct dependencies are pinned in `pyproject.toml`: Astropy 8.0.2, NumPy 2.5.4, Matplotlib 3.11.2, JupyterLab 4.6.6, ipykernel 7.4.0, nbconvert 7.17.2; pytest 9.2.0 is used for checks. `uv.lock` also pins transitive dependencies, including ERFA and time tables. Installation requires access to packages; subsequent calculations work offline without catalogue requests. The script disables automatic IERS downloads.

One date, ISO in UTC by default:

```sh
uv run --offline --locked python barnard_motion.py \
  --date 2025-04-01T00:00:00 \
  --output output/barnard_single.csv
```

Numeric dates require an explicit format **and** time scale:

```sh
uv run --offline --locked python barnard_motion.py \
  --date 60800.0 --format mjd --scale utc \
  --output output/barnard_mjd.csv

uv run --offline --locked python barnard_motion.py \
  --date 2016.0 --format jyear --scale tcb \
  --output output/barnard_reference_epoch.csv
```

The last example checks the original epoch. MJD 60800.0 above is a separate input example; it **does not represent the same instant** as April 1.

Supported formats: `iso`, `isot`, `mjd`, `jyear`; scales: `utc`, `tai`, `tt`, `tcg`, `tdb`, `tcb`. ISO accepts `YYYY-MM-DD`, `YYYY-MM-DDTHH:MM:SS[.fraction]`, and a space in place of T. A calendar date alone means 00:00:00 in the selected scale. The Z suffix is allowed only for UTC. The computer time zone is not used; entries with time-zone offsets are rejected. `jyear` is a Julian year: 2000.0 corresponds to JD 2451545.0 in the selected scale and a year is 365.25 days; it is not a decimal calendar year. UT1 is not supported in this calculation.

Full list of arguments:

```sh
uv run --offline --locked python barnard_motion.py --help
```

Without `--output`, the result is written to `output/barnard_positions.csv`. `--catalog` explicitly specifies the parameter CSV path. By default the catalogue is sought beside the script, regardless of the current directory. An output file of the same name is overwritten. On invalid input the program exits with code 2 and an explanation.

## Observation list

Tutorial run:

```sh
uv run --offline --locked python barnard_motion.py \
  --observations data/barnard_observations_demo.csv \
  --output output/barnard_demo_positions.csv
```

The example has exactly two dates: **2025-04-01 00:00 UTC** and **2025-10-01 00:00 UTC**. They are marked `is_demo=true`; `obs_id`, detector, release and link are empty. These are not discovered frames. For real data, fill in a file using the header in `data/barnard_observations_template.csv`.

| Input column | Contents |
| --- | --- |
| `obs_id` | Actual observation ID as a string. A tutorial row may be empty when `is_demo=true`. |
| `time_format`, `time_scale` | Required in every row, including ISO. Not inferred from the date value or filename. |
| `time` | A single instant in the stated format and scale. |
| `t_min`, `t_max` | Instead of `time`, the interval start and end in the same format and scale. |
| `detector`, `data_release` | Optional detector number/name and release; preserved unchanged. |
| `access_url` | Optional frame or DataLink resource link; retained but not opened. |
| `is_demo` | `true` or `false`, default `false`. |

Different rows in one CSV may use different formats and scales. Providing both `time` and bounds, incomplete bounds, a reversed interval, or no time causes rejection with the row number. Equal bounds are allowed. The midpoint is calculated as `t0.tai + (t1.tai - t0.tai)/2`, accounting for UTC leap seconds rather than averaging strings or calendar hours. It is the unweighted interval midpoint, not the effective time of a flux-weighted measurement.

`--format` and `--scale` apply only to `--date`; in batch mode they are read from the CSV. Batch mode expects normalized CSV, not IPAC `.tbl`. When transferring from IRSA, explicitly set the verified time scale, transfer `obs_collection` to `data_release`, and the detector number to `detector`. `obs_id` may repeat for different detectors; rows are not combined.

Identifiers are read as strings by a standard CSV parser, preserving leading zeros and all Gaia digits. When opening the output in Excel or Sheets, import `obs_id` and `source_id` as **text** to prevent automatic rounding of long numbers.

## Tutorial notebook

```sh
uv run --offline --locked jupyter lab output/jupyter-notebook/barnard_motion.ipynb
```

Use the Python kernel from `.venv`. The notebook shows the original parameters and correlations, the selected date and result, a comparison of two dates, a plot, and an exercise changing the date. To reproduce everything, restart the kernel and run all cells. The notebook imports `barnard_motion.py`; it does not contain a second propagation implementation.

Run the entire notebook without the interface:

```sh
uv run --offline --locked jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=180 output/jupyter-notebook/barnard_motion.ipynb
```

The notebook saves the same `output/barnard_demo_positions.csv` as the script, plus `output/barnard_motion_demo.png`. The executed notebook and tutorial results are saved in the repository. Rerunning it updates both output files.

## Method and source record

Only the exact string `source_id=4472832130942575872` is selected from the CSV, checking Gaia DR3 release, ICRS, epoch J2016.0 TCB, the RA error convention, and required values. Mixing rows or substituting zero for a missing RV is not permitted.

| Parameter | Value | Formal 1σ error |
| --- | ---: | ---: |
| RA at J2016.0 TCB, degrees | 269.44850252543836 | σα⋆ = 0.026238997 mas |
| Dec at J2016.0 TCB, degrees | 4.739420051112412 | 0.029043527 mas |
| μα⋆, mas/year | −801.5509783684709 | 0.031820867 |
| μδ, mas/year | 10362.394206546573 | 0.036070455 |
| Parallax, mas | 546.975939730948 | 0.040116355 |
| RV, km/s | −110.46822 | 0.13125522 |

The official [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html) defines the epoch in TCB, RA/Dec as barycentric ICRS, and `pmra` and `ra_error` with the cos δ factor. The [Gaia MAST page](https://archive.stsci.edu/missions-and-data/gaia) is source resource #4; the selected record and reproducible ESA query are already in the CSV. Links to the record are also included in each result.

`SkyCoord` receives `pm_ra_cosdec`, `pm_dec`, distance from parallax, RV, and the reference epoch. `SkyCoord.apply_space_motion(new_obstime=time.tcb)` propagates the position using a constant spatial velocity model. Unlike simply adding μα⋆ and μδ, the changing distance causes a perspective change in angular motion. Implementation and model conventions are defined by the [Astropy method](https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html) and the ERFA version used.

For uncertainties, the local perturbation vector **q = (δα⋆, δδ, δϖ, δμα⋆, δμδ, δRV)** is used, with units mas, mas, mas, mas/year, mas/year, km/s. The five-parameter covariance is reconstructed from the ten Gaia correlations as **Cᵢⱼ = ρᵢⱼ σᵢ σⱼ**. Independent RV is added: C₆₆ = σ²RV, with all other elements of its row/column zero. Values are checked for finiteness, errors for positivity, correlations for valid bounds, and the matrix for positive semidefiniteness.

A numerical 2×6 Jacobian J is calculated by central differences of the same propagation function. Input perturbation δα⋆ is converted to an RA change by dividing by cos δ at the **original** position. Output offsets use a local spherical east/north frame relative to the central propagated position; its differentials are dα cos δ and dδ **at the output epoch**. **C_out = J C_in Jᵀ** is a 2×2 position-only covariance, not a complete new six-parameter Gaia solution.

Difference steps: 10 mas for each coordinate, 0.1 mas for parallax, 0.1 mas/year for each motion component, and 0.1 km/s for RV. These exceed angular machine-rounding noise and are small enough for a local derivative. Checks with half and double steps test output covariance stability. `apply_space_motion` does not automatically propagate uncertainties.

## Output CSV

| Columns | Meaning |
| --- | --- |
| `obs_id`, `detector`, `data_release`, `access_url`, `is_demo` | Link to input row; `obs_id` is empty for a single date. |
| `input_time_format`, `input_time_scale`, `input_time`, `input_t_min`, `input_t_max` | Format/scale and original values for reproduction. |
| `time_method` | `exact` or `midpoint_tai`. |
| `time_utc`, `time_tcb` | Instant used, in ISO with up to 9 fractional-second digits. Z is used only for UTC; TCB is indicated by the column name. |
| `time_mjd_utc`, `time_jyear_tcb` | Additional numeric representations of the same instant. |
| `coordinate_frame`, `position_type` | `ICRS`, `barycentric_astrometric`. |
| `ra_deg`, `dec_deg` | Propagated coordinates in degrees, without annual parallax. |
| `sigma_alpha_star_mas`, `sigma_delta_mas` | Formal 1σ standard errors: √C₀₀, √C₁₁. The RA error includes cos δ. |
| `ra_dec_corr` | C₀₁ / (σα⋆ σδ), correlation at the output epoch. |
| `cov_alpha_star_delta_mas2` | C₀₁ in mas². Together with the errors, permits reconstruction of the 2×2 matrix. |
| `distance_pc` | Propagated distance in parsecs. |
| `parallax_scale_1au_mas` | 1000 / `distance_pc`, characteristic parallax scale for a 1 AU baseline, in mas. |
| `parallax_applied` | Always `false`; this scale is not added to RA/Dec. |
| `source_id`, `gaia_source_url` | Complete exact ID and link to the original Gaia record. |

Rounded tutorial results, for reference:

| Date, 00:00 UTC | RA, degrees | Dec, degrees | σα⋆, mas | σδ, mas | 1 AU scale, mas |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2025-04-01 | 269.446435307 | 4.766052499 | 0.290077 | 0.324998 | 547.288714 |
| 2025-10-01 | 269.446323228 | 4.767496369 | 0.305989 | 0.344398 | 547.305668 |

Over the interval the local displacement relative to the first position is **−0.402086″ east** (i.e. west) and **+5.197935″ north**. These are spherical displacement components, practically equal here to Δα cos δ at the first point and Δδ. Formal errors at different dates are correlated through the common source record; the error of the difference is not the square root of the sum of independent variances. The script outputs errors for each position, not the error of their difference.

## Provided files

The supplied files `/Users/dmitrytukhvatullin/Desktop/table_Spectral-Images (1).tbl` and `/Users/dmitrytukhvatullin/Desktop/table.tbl` were checked. Their contents are used as data, not instructions.

- The first file contains 372 rows and 188 distinct `obs_id` values; collections `spherex_dr1` (265 rows) and `spherex_qr3` (107), and detectors D1–D6. Different detectors for the same observation have separate rows.
- Frame centres: RA **294.32186906–298.26898590°**, Dec **−16.67792059…−12.85175631°**. All coverage polygon vertices are within RA 292.3245–300.3479°, Dec −18.6877…−10.7726°. Barnard's Star is near RA 269.446°, Dec +4.77° on the dates in question; the nearest centre is approximately **30.68°** away. These frames do not cover Barnard's Star.
- `table.tbl` is a metadata record, not another frame list: 62 fields for **2025W38_1A_0596_1**, **SPHEREx-D6**, **spherex_dr1**. Its coordinates and interval match the first row of the larger file.
- The record defines `t_min`/`t_max` as **MJD**, but **the supplied files do not specify a time scale**. The complete range is MJD 60935.405692–61317.775151; the scale must be confirmed from product documentation or a FITS header before calculation. Product creation and release times are not observation times.
- Every `access_url` in the larger file has MIME type `application/x-votable+xml;content=datalink`: these are **DataLink** links, not direct FITS links. They help locate products, but the FITS files themselves were not downloaded here.

The rows in these files were therefore **not used as real observations of Barnard's Star**. The next application requires a team list for its actual region: ID, FITS/cutout link or resolvable DataLink, start and end with format and scale, detector, and release. The search structure is shown in the [official IRSA/SPHEREx introduction](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/). One suitable frame suffices to check the structure; comparing motion requires different dates.

## Verification

```sh
uv run --offline --locked pytest -q
```

The tests check:

- the reference epoch: coordinate difference <0.001 mas, errors and covariance within 1% of the catalogue;
- direction and magnitude of motion over a short interval, including a separate cos δ check at declination 60°;
- agreement of short-term covariance with analytical position/motion propagation, and use of correlations;
- the influence of RV on position and of its independent uncertainty on output covariance;
- numerical Jacobian stability when step sizes change;
- equivalence of ISO UTC, MJD UTC, and Julian year/ISO TCB representations of the same instant;
- the midpoint across a UTC leap second; validation of ambiguous or reversed intervals, invalid dates, scales, and parameters;
- agreement of single, batch, and direct function calls, and exact preservation of long IDs and metadata;
- calculation with network connections disabled.

The notebook was also run completely from a fresh kernel and checks tutorial results against the shared functions. Completion is checked by the absence of errors in every cell; the tutorial CSV matches a separate script run.

## Scientific limitations

1. The constant spatial velocity model includes perspective changes in angular motion, but not arbitrary physical accelerations or orbital motion. A small formal error does not justify distant extrapolation.
2. Errors are a local linear approximation of formal Gaia errors and independent RV error. Gaia/RV systematics, parallax zero point, time uncertainty, spatial registration errors, and centroid measurement errors are excluded. The source record's RUWE is ≈1.08485, but this does not turn formal errors into a complete accuracy estimate.
3. Distance is obtained by simply inverting the raw positive parallax; no Bayesian distance inference or zero-point correction is applied. The relative parallax error is small for this record. The 1 AU scale characterises distance; it is not a prediction of an individual parallactic displacement or its uncertainty.
4. The SPHEREx observer position, annual/topocentric parallax, aberration, and gravitational deflection for an image are absent. The prediction cannot be treated as a ready apparent detector position. Comparison with a particular image requires the observer position and checks of its WCS conventions.
5. The supplied instant is used without an additional barycentric photon-arrival-time correction for a particular observer. Such a correction does not follow from UTC↔TCB conversion alone. Time and geometry must be checked separately for precise comparison with a frame.
6. Time-scale conversion uses pinned ERFA and leap-second tables. UTC for distant future dates depends on leap seconds that have not yet been announced; ERFA dubious-year warnings are not suppressed. Changing the computer's time zone does not affect the calculation.

IRSA search, FITS download and reading, quality flags, precise observer treatment, and comparison with measured coordinates remain the team's next stage. WiseView and NEOSSat/CADC are not needed for this local calculation.
