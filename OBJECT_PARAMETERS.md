# Verified astrometric parameters for five targets

Checked on **November 14, 2026**, Almaty time. The numbers were checked against Gaia Archive records and the original tables in the publications. Verification against a source does not mean that the model accounts for every feature of an object.

Machine-readable data: [OBJECT_PARAMETERS.csv](OBJECT_PARAMETERS.csv). It contains exactly five rows, the original numerical precision, and separate uncertainty columns. Gaia values below are rounded for readability. All coordinates are in **ICRS** at the stated epoch; positions at the dates of the SPHEREx frames have not been calculated here.

## Coordinates and epochs

| Target; what the record describes | RA, degrees | Dec, degrees | Position epoch | RA error, mas | Dec error, mas |
| --- | ---: | ---: | --- | ---: | ---: |
| Barnard's Star | 269.448502525 | +4.739420051 | J2016.0, TCB | 0.026239 | 0.029044 |
| Epsilon Indi B: Gaia source identified with **Ba** | 331.076452597 | −56.793812067 | J2016.0, TCB | 0.504801 | 0.467587 |
| WISE 0855−0714 | 133.780984 | −7.243932 | **MJD 57633.39**, time scale not specified in the table | 2.3¹ | 2.3 |
| RW Aur: **A as identified in the literature**, eastern component | 76.956544311 | +30.401316150 | J2016.0, TCB | 0.281577 | 0.177497 |
| R Dra | 248.167378367 | +66.755006745 | J2016.0, TCB | 0.066872 | 0.066013 |

For Gaia rows, the RA error is **σα⋆ = σα cos δ**, not the error in the angle α without the cos δ factor. For WISE, the published RA error is retained; note 1 explains the distinction. The epoch and errors are defined in the [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html).

## Proper motion and parallax

| Target | μα⋆ = (dα/dt) cos δ, mas/year | μδ, mas/year | Parallax, mas | Gaia RUWE |
| --- | ---: | ---: | ---: | ---: |
| Barnard's Star | −801.550978 ± 0.031821 | +10362.394207 ± 0.036070 | 546.975940 ± 0.040116 | 1.084851 |
| Epsilon Indi B: Gaia record for Ba | +3981.976663 ± 0.601657 | −2466.831815 ± 0.625683 | 270.658032 ± 0.689557 | **4.295541** |
| WISE 0855−0714 | −8123.7 ± 1.3 | +673.2 ± 1.3 | 439.0 ± 2.4 | — |
| RW Aur A as identified in the literature | +4.906034 ± 0.329158 | −24.577535 ± 0.254649 | 5.454096 ± 0.332450 | **16.422132** |
| R Dra | −20.642812 ± 0.085128 | +6.148819 ± 0.091321 | 1.250918 ± 0.064788 | 1.297982 |

All five astrometric parameters in each row come from one solution. For the four Gaia rows, the source is `gaiadr3.gaia_source`; for WISE, it is Table 5 of Kirkpatrick et al. (2021). Gaia errors are formal standard errors, 1σ. For WISE, the published fit errors are retained: Table 5 does not separately specify a confidence level. The errors have not been scaled, and no Gaia parallax zero-point correction has been applied.

The CSV also contains Barnard's Star's radial velocity of **−110.46822 ± 0.13125522 km/s** from Gaia DR3. The velocity fields for the other four rows are empty because the value is absent from the selected record or table. An empty field does not mean zero velocity or a lack of measurements in other works.

## Identifiers and sources

| Target | Exact record used for the numbers | Identity check |
| --- | --- | --- |
| Barnard's Star | [Gaia DR3 **4472832130942575872**](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=I/355/gaiadr3&Source=4472832130942575872) | [SIMBAD: NAME Barnard's star](https://simbad.harvard.edu/simbad/sim-id?Ident=Barnard%27s+star); GJ 699, HIP 87937 |
| Epsilon Indi B | [Gaia DR3 **6412596012146801152**](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=I/355/gaiadr3&Source=6412596012146801152) | [SIMBAD: * eps Ind B](https://simbad.harvard.edu/simbad/sim-id?Ident=%2A+eps+Ind+B); pair WDS J22034−5647Ba,Bb; Gaia source assigned to Ba by [Baig et al. (2024), Table 2, p. 3790](https://uhra.herts.ac.uk/id/eprint/11346/1/stae2005.pdf#page=7) |
| WISE 0855−0714 | [Kirkpatrick et al. (2021), Table 5](https://arxiv.org/html/2011.11616#S4.T5); [VizieR J/ApJS/253/7/plx, PName=0855−0714](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJS/253/7/plx&PName=0855-0714) | [SIMBAD: WISE J085510.83−071442.5](https://simbad.harvard.edu/simbad/sim-id?Ident=WISE+J085510.83-071442.5); SIMBAD's main identifier is WISEA J085510.74−071442.5 |
| RW Aur A as identified in the literature | [Gaia DR3 **156430817820015232**](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=I/355/gaiadr3&Source=156430817820015232) | [Current SIMBAD entry: V* RW Aur B](https://simbad.harvard.edu/simbad/sim-id?Ident=V%2A+RW+Aur+B); the component-label ambiguity is discussed below |
| R Dra | [Gaia DR3 **1648328016354496384**](https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=I/355/gaiadr3&Source=1648328016354496384) | [SIMBAD: V* R Dra](https://simbad.harvard.edu/simbad/sim-id?Ident=R+Dra); HD 149880, HIP 81014 |

The four Gaia records were retrieved directly from [ESA Gaia Archive TAP](https://gea.esac.esa.int/tap-server/tap) by querying `gaiadr3.gaia_source` using integer `source_id` values. In the CSV, `source_url` contains a reproducible query for each record; the VizieR links above provide a convenient view of the catalogue. Release description: [Gaia Collaboration et al. (2023)](https://doi.org/10.1051/0004-6361/202243940); astrometric solution: [Lindegren et al. (2021)](https://doi.org/10.1051/0004-6361/202039709).

The WISE row was read from the [machine-readable VizieR table](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J%2FApJS%2F253%2F7%2Fplx&PName=0855-0714&-out.all=&-out.max=10) and independently checked against Table 5 of [Kirkpatrick et al. (2021)](https://doi.org/10.3847/1538-4365/abd107). The fields used were `RA_ICRS`, `DE_ICRS`, `MJD`, `plx`, `pmRA`, `pmDE`, and their corresponding errors. The VizieR-computed columns `_RA.icrs` and `_DE.icrs` at epoch J2000 were not used.

## Object-specific limitations

### Epsilon Indi B: component and barycentre

The team's target is the **Ba/Bb** pair, while the main CSV row describes the Gaia source that [Baig et al. (2024)](https://doi.org/10.1093/mnras/stae2005) associate with **Ba**. This is not a measurement of the Ba/Bb barycentre, nor is it a guaranteed position of the pair's photocentre in the SPHEREx band. Bb's possible influence on the optical photocentre cannot be converted into a precise correction without an additional model. The elevated RUWE also calls for caution.

A separate system solution from [Chen et al. (2022), Table 4 and §6.2](https://arxiv.org/html/2205.08077#S6.T4), which accounts for the Ba/Bb orbit, gives:

| Parameter | Chen et al. (2022) system solution, posterior mean ± 1σ |
| --- | ---: |
| μα⋆, mas/year | +3987.41 ± 0.12 |
| μδ, mas/year | −2505.35 ± 0.10 |
| Parallax, mas | 274.99 ± 0.43 |

These values have **not been substituted** alongside the Gaia coordinates: they describe a different model. Table 4 does not contain an absolute barycentre RA/Dec pair with an epoch and errors sufficient to replace the complete row. In subsequent calculations, use the Gaia row only with its limitation explicitly stated; precise system motion will require a consistent system solution. Epsilon Indi A parameters are not used here.

### RW Aur: A/B ambiguity

**In this table, A means the eastern component identified as A in studies of the resolved pair.** This choice is based on coordinates and geometry, not solely on the letter in a catalogue entry.

On the date checked, [SIMBAD labels Gaia DR3 156430817820015232 as B](https://simbad.harvard.edu/simbad/sim-id?Ident=V%2A+RW+Aur+B), and [156430822114424576 as A](https://simbad.harvard.edu/simbad/sim-id?Ident=V%2A+RW+Aur+A). However, [Dodin et al. (2019), note 1](https://doi.org/10.1093/mnras/sty2988) explicitly assign the second identifier to **B**. [Kurtovic et al. (2024), §3.1](https://arxiv.org/html/2407.18828v1#S3.SS1) place B southwest of A. In the two Gaia records at J2016.0, the second source is indeed 1.438″ west and 0.402″ south of the first: separation ≈1.493″, position angle ≈254.40°. These angles were calculated during the check; they are not a new observation.

Therefore, **156430817820015232** was selected for the main row. Its RUWE of 16.422132 is consistent with the discussion in [Kurtovic et al. (2024), §4.6](https://arxiv.org/html/2407.18828v1#S4.SS6) of a RUWE of about 16 for A, binarity, variability, and Gaia systematics. The formal errors in the row do not cover all these effects. The parameters of the southwestern B component do not replace those of A; when the pair is unresolved, the measured centre and flux may depend on both components.

### WISE 0855−0714

All parameters refer to a single five-parameter fit in Table 5 at **MJD 57633.39**. This is the reference position of the motion model, with parallax modelled separately, not an observed position from a particular frame. The fit used 19 Spitzer epochs and 46 Gaia DR2 reference stars; the row has quality category `h`, χ²=19.56 for 33 degrees of freedom, and χ²/ν=0.59. See [Table 5 and the methodology of Kirkpatrick et al. (2021)](https://arxiv.org/html/2011.11616#S4.T5).

¹ Table 5 calls both position errors uncertainties in mas, but does not explicitly specify the projection of the RA error. Therefore, `ra_error_basis` differs from Gaia's; the value 2.3 has not been converted or automatically designated σα⋆. The MJD time scale has likewise not been assumed. The coordinates were published to six decimal places; no additional significant digits have been reconstructed. Radial velocity and parameter correlations are absent from the selected table.

The RA proper motion already includes cos δ, as follows from equation (3) in the methodology of [Kirkpatrick et al. (2019), §5.2.3](https://arxiv.org/html/1812.01208#S5.SS2.SSS3), referenced by the 2021 paper. There is no need to divide or multiply this value by cos δ when reading the CSV.

### Barnard's Star and R Dra

Barnard's Star has a large proper motion and an available RV; accurate propagation requires a spatial motion model, including perspective acceleration. The RV in the CSV is taken from Gaia, separately from the five astrometric parameters. R Dra's identity as a Mira variable has been confirmed; catalogue errors alone do not assess all possible photocentre shifts of a variable star. [SIMBAD: Barnard's Star](https://simbad.harvard.edu/simbad/sim-id?Ident=Barnard%27s+star), [R Dra](https://simbad.harvard.edu/simbad/sim-id?Ident=R+Dra).

## How to read the CSV

Format: UTF-8, comma-separated, decimal point, with column headers in the first row. **Read `source_id` as text**: automatically converting long Gaia identifiers to floating-point numbers may corrupt their final digits. Empty cells indicate a missing value in the selected source or a field that does not apply, not zero.

| Fields | Meaning |
| --- | --- |
| `object`, `canonical_id`, `component`, `astrometric_target` | Target, exact identifier, component, and physical meaning of the main row |
| `simbad_id`, `identity_url` | Current SIMBAD name and identity-check link; for RW Aur, the component letter differs from the literature |
| `coordinate_frame` | ICRS; this denotes the reference system, not the position epoch |
| `ref_epoch`, `ref_epoch_format`, `ref_epoch_scale` | Epoch value, `jyear`/`mjd` format, and time scale. For Gaia: 2016.0, `jyear`, TCB. For WISE: 57633.39, `mjd`, time scale empty |
| `ra_deg`, `dec_deg` | Original coordinates at the reference epoch, in degrees |
| `ra_error_mas`, `ra_error_basis`, `dec_error_mas` | Position errors and definition of the RA error; see note 1 |
| `pmra_mas_yr`, `pmdec_mas_yr` and corresponding `*_error_*` fields | μα⋆ and μδ with errors in mas/year; positive signs indicate increasing RA and Dec |
| `parallax_mas`, `parallax_error_mas` | Published absolute parallax and its error in mas |
| `radial_velocity_km_s`, `radial_velocity_error_km_s` | RV from the selected Gaia record, if available; populated only for Barnard's Star here |
| `catalog`, `source_id`, `source_bibcode`, `source_url`, `checked_date` | Origin of the numbers, record identifier, bibliographic code, link, and date checked |
| `ruwe`, `astrometric_params_solved`, `astrometric_excess_noise_mas` | Gaia solution diagnostics; 31 means five parameters, 95 means six including pseudocolour, not an orbital solution |
| `fit_chi2`, `fit_dof`, `fit_reduced_chi2` | WISE fit diagnostics; these values are not equivalent to RUWE |
| Ten `*_corr` fields | Original correlations between the five Gaia astrometric parameters; empty for WISE |
| `uncertainty_description`, `notes` | Definition of uncertainties, missing fields, and limitations on use |

Gaia correlations are retained for future uncertainty propagation: the parameter order is **(α⋆, δ, ϖ, μα⋆, μδ)**; covariance equals `corr(i,j) × sigma(i) × sigma(j)`. This is the submatrix for the five astrometric parameters, not the full matrix including pseudocolour for type 95 solutions. Formal errors and correlations do not replace an assessment of systematic effects in binary and variable objects. [Gaia DR3 field definitions](https://irsa.ipac.caltech.edu/data/Gaia/dr3/gaia_dr3_source_colDescriptions.html).

## Checks performed

- Checked SIMBAD and Gaia identifiers; investigated the RW Aur A/B ambiguity using independent publications and relative positions.
- Transferred Gaia parameters from the original response without rounding or combining them with other solutions; checked the WISE row against VizieR and Table 5.
- Preserved the actual reference epoch for all coordinates. SIMBAD coordinates displayed at J2000 were not mixed with the original Gaia coordinates at J2016.0.
- Checked units and directions of motion: Barnard's Star moves mostly north, Epsilon Indi Ba southeast, and WISE west and slightly north.
- Re-read the CSV: five unique targets, valid numeric fields, complete identifiers, and no missing values replaced with zeros.

Image retrieval, propagation to observation dates, photometry, and determination of the methane index are outside the scope of these files.
