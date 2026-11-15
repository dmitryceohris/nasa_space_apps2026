# Scientific acceptance of targets and selected observations

Date: **November 14, 2026, Almaty**. Base: `origin/main` at **`eaca1e5752f9d7faf2afd94b49df131c91cadab1`**, checked after `git fetch`. Working branch: `codex/science-acceptance`.

**Decision:** target identities and display of the selected fields are accepted within the limits below. The saved images do not establish detection of all eight sources, or measurements of proper motion, variability, or period. This review does not validate the scientific fitness of the centroids or photometry.

## Evidence and scope

Reviewed [OBJECT_PARAMETERS.md](../../OBJECT_PARAMETERS.md), [OBJECT_PARAMETERS.csv](../../OBJECT_PARAMETERS.csv), the since-deleted input-data assessment, [TARGET_VALIDATION.md](../TARGET_VALIDATION.md), and [DILAYLA_ACCEPTANCE.md](../DILAYLA_ACCEPTANCE.md). The main sample is [public-final](../../data/validation/public-final/target_cases.csv), not the intermediate public-preview: **43 frames, 43 distinct obs_ids, one product per frame**, D6 for WISE and D3 for the others, `single`, FOV=6′. Four frames from the [exact WISE coordinate query](../../data/validation/responses/wise-objective.json) were checked separately. Browser evidence in [browser-public-final.json](../../data/validation/browser-public-final.json) is accepted as a check of image loading and date switching; the browser was not rerun.

New evidence is in [data/validation/science-acceptance](../../data/validation/science-acceptance/):

- [sources.json](../../data/validation/science-acceptance/sources.json) — URLs, retrieval dates, SHA-256 hashes of primary responses, and retrieval limitations.
- [gaia.csv](../../data/validation/science-acceptance/gaia.csv) — raw ESA TAP response, including the four rows used and the neighbouring RW Aur record. The values, errors, diagnostics, and ten correlations in the four rows **exactly matched** the existing OBJECT_PARAMETERS.csv. This preserves the source response; it does not rebuild the parameter table.
- [wise-plx.tsv](../../data/validation/science-acceptance/wise-plx.tsv) — VizieR response for the original Kirkpatrick table; the epoch and parameters matched the CSV.
- [rw-A-simbad.txt](../../data/validation/science-acceptance/rw-A-simbad.txt), [rw-B-simbad.txt](../../data/validation/science-acceptance/rw-B-simbad.txt), and [rdra-simbad.txt](../../data/validation/science-acceptance/rdra-simbad.txt) — new SIMBAD responses checking component letters, identifiers, and variable-star type.
- [irsa-products.csv](../../data/validation/science-acceptance/irsa-products.csv) — 43 selected URIs from the existing IRSA metadata cache. Full product paths, including `dr1`/`qr3` and processing version, are preserved. Fresh IRSA TAP requests timed out; **this is not a new archive response or a completeness check**.
- [evidence.json](../../data/validation/science-acceptance/evidence.json) — mapping of each response and product, SHA-256 hashes of responses and available FITS cutouts, original header times, WCS position, FLAGS values in the 5×5 region, and pixel counts for each set bit. Full FITS files were not added to Git.

Matching used **exact URI + OBSID + DETECTOR**, not date or a similar name. All 47 checked cutouts have `RADESYS=ICRS` and `TIMESYS=UTC`; time comes from `MJD-AVG` and the UTC date matched the API. Rounding MJD to four decimals in the API produces a maximum rounding error of 4.32 s; exact values remain in evidence.json. The lower bound of the archive interval was not presented as the exposure midpoint.

**34 DR1 products and 9 QR3 products** belong to public-final; the four WISE objective products are the same four DR1 products with different cutouts. IRSA describes DR1 and QR3 as a jointly processed set; this permits viewing them together but does not guarantee identical bands or the fitness of every pixel. [IRSA DR1 Overview](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_dr1.html), [IRSA Data Access](https://irsa.ipac.caltech.edu/spherex-archive-documentation/spherex-data-access/).

## Verdict for each target

“Accepted” here describes this scientific review; it does not claim agreement from the owners of shared files. A detailed list of all obs_ids, releases, dates, bands, and flags follows.

| Target | Source / component | Checked observations | Quality / limitations | Supported conclusion | Verdict |
| --- | --- | --- | --- | --- | --- |
| Barnard’s Star | Gaia DR3 **4472832130942575872**, GJ 699 / HIP 87937; J2016.0 TCB | [Response](../../data/validation/public-final/responses/0-barnard.json): 5 dates, 2025-08-15—2026-08-16; D3, 4 DR1 + 1 QR3 | All centroids `poor`; OVERFLOW, SUR_ERROR, TRANSIENT and BLOOM near the prediction on every date; λ=1.6504–2.3074 μm | Field coverage and expected northwest motion; numerical displacement is a prediction | **Accepted for viewing and prediction; motion measurement remains open** |
| Epsilon Indi B | Gaia DR3 **6412596012146801152**, literature **Ba** of the Ba/Bb system; not Epsilon Indi A | [Response](../../data/validation/public-final/responses/1-epsindib.json): 7 dates, 2025-04-29—2026-09-22; D3, 6 DR1 + 1 QR3 | RUWE=4.2955413; Ba/Bb unresolved; all 7 results `ambiguous`. NONFUNC near prediction on first date, only SOURCE on the others | The B-system field and Ba reference may be shown with the model caveat; the reference is not the barycentre or IR photocentre of the pair | **Component accepted; precise system astrometry requires an orbit and flux model** |
| WISE 0855−0714 | WISE J085510.83−071442.5; one Kirkpatrick Table 5 row, **MJD 57633.39** | [Preset](../../data/validation/public-final/responses/2-wise0855.json) and [exact coordinates](../../data/validation/responses/wise-objective.json): 4 dates each, 2025-05-05—2026-04-21; identical obs_ids, D6 / DR1 | Preset λ=4.4178–4.9674 μm; 3 `ambiguous`, 1 `ok/fair`; exact query has 4 `ambiguous`, old-schema quality fields absent; detection not confirmed | Both fields are available; the model lies within each 6′ field, including the fixed objective centre | **Coverage accepted; detection and motion require additional data** |
| RW Aur | Gaia DR3 **156430817820015232**, eastern **A in the literature**; SIMBAD currently labels this record B | [Response](../../data/validation/public-final/responses/3-rwaur.json): 3 dates, 2025-09-11—2026-09-13; D3, 2 DR1 + 1 QR3 | A/B separation ≈1.493″ versus 6.15″ pixel; RUWE=16.422132. `poor`, `fair`, `fair`; masks near prediction on all three dates, OVERFLOW on two | View the unresolved system with the literature A position as a reference; differences in bands/fluxes do not establish a change in A itself | **Identity accepted with an explicit caveat; separate A photometry not accepted** |
| R Dra | Gaia DR3 **1648328016354496384**, V* R Dra / HD 149880 / HIP 81014; Mira | [Response](../../data/validation/public-final/responses/4-rdra.json): 11 distinct dates, 2025-04-24—2026-06-12; D3, 9 DR1 + 2 QR3 | All 11 centroids `poor`; OVERFLOW/SUR_ERROR/BLOOM on all dates; λ=1.6297–2.4032 μm | A real sequence of different dates exists; brightness comparisons require suitable comparable bands | **Viewing accepted; measured variability and period not confirmed** |
| Proxima Cen | Sesame/SIMBAD **NAME Proxima Centauri**, GJ 551 / HIP 70890; saved aliases include Gaia DR3 5853498713190525696 | [Response](../../data/validation/public-final/responses/5-proxima-cen.json): 4 dates, 2025-07-22—2026-07-20; D3, 3 DR1 + 1 QR3; name and coordinates return identical frames | Fixed catalogue position, epoch unverified; all `ambiguous`; PERSIST in one neighbourhood of the fixed centre | Resolver, field and sequence checked; matching queries do not establish a measurement of Proxima's motion | **Accepted for a fixed field** |
| 61 Cyg A | Sesame/SIMBAD **\* 61 Cyg A**, HD 201091 / HIP 104214 / GJ 820 A; Gaia DR3 alias 1872046609345556480 | [Response](../../data/validation/public-final/responses/6-61-cyg-a.json): 5 dates, 2025-05-25—2026-07-05; D3, 3 DR1 + 2 QR3; name and coordinates equivalent | Epoch unverified, position fixed; all `ambiguous`; λ=1.6831–2.2257 μm; possible contribution from neighbouring B | The canonical name refers to A; this does not confirm the component associated with every detected spot | **Identifier and field accepted; candidate association remains open** |
| Vega | Sesame/SIMBAD **\* alf Lyr**, NAME Vega / HD 172167 / HIP 91262 | [Response](../../data/validation/public-final/responses/7-vega.json): 4 dates, 2025-08-28—2026-08-30; D3, 3 DR1 + 1 QR3; name and coordinates equivalent | Epoch unverified; all centroids `poor`; OVERFLOW/SUR_ERROR in every fixed-centre neighbourhood, BLOOM | View the field of a bright source and demonstrate mask limitations | **Identifier and field accepted; precise astrometry/photometry not accepted** |

## Physical interpretation

### Barnard’s Star

Negative μα⋆ and positive μδ mean motion west and north. For the first and last public-final dates, recalculation with the existing `sky_motion` gives geocentric offsets **east=−0.80794″, north=+10.38820″** over 365.7661 days. On a north-up grid with east to the left, this is movement right and up: saved PNG markers `[28.855, 31.044] → [28.986, 29.355]` have the correct signs. This is a catalogue calculation, not measured centroids. RV=−110.46822 km/s contributes to spatial propagation and perspective change in motion. [Gaia DR3 definitions](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html).

The original supplied export concerns the NGC 6822 field and does not replace GJ 699 observations. Its initial D1/DR1 series `2025W35_1A_0432_4`, `2025W40_1A_0471_1`, `2026W11_1A_0539_2` and barycentric prediction use **a different sample and direction model**. Its 5.594198″ must not be presented as measured displacement or as the result of the public-final D3 series. Matching detector numbers also do not imply identical local bands.

### Epsilon Indi B: Ba, barycentre and photocentre

[Baig et al. (2024), Table 2, p. 3790](https://uhra.herts.ac.uk/id/eprint/11346/1/stae2005.pdf#page=7) associate 6412596012146801152 with **eps Indi Ba**, not A or Bb alone. [Chen et al. (2022), Table 4](https://arxiv.org/html/2205.08077#S6.T4) give a system solution accounting for orbit: μα⋆=3987.41, μδ=−2505.35 mas/yr, ϖ=274.99 mas; it cannot be substituted alongside the Gaia Ba position without a consistent reference position and epoch. This is a difference between solutions, not CSV rounding error.

The barycentre is mass-weighted, while the photocentre is weighted by flux in a particular band: `r_bar = (mBa*rBa + mBb*rBb)/(mBa+mBb)`; `r_ph(λ) = (FBa(λ)*rBa + FBb(λ)*rBb)/(FBa(λ)+FBb(λ))`. Thus Ba's position, the barycentre, and the Ba/Bb photocentre differ. Without the orbit and flux ratio as a function of λ, there is no exact correction. The single Gaia record's prediction is a useful reference; it does not model subsequent orbital motion. SOURCE at the expected position does not establish detection of Ba; `good` for one ambiguous candidate likewise does not resolve its identity.

### WISE: two independent modes and coverage check

In [Kirkpatrick et al. (2021), Table 5](https://arxiv.org/html/2011.11616#S4.T5), coordinates **133.780984°, −7.243932°** refer to **MJD 57633.39**, not J2000 or the SPHEREx epoch. Proper motion is −8123.7, +673.2 mas/yr and parallax 439.0 mas; parameters and published errors come from one fit (19 Spitzer epochs, 46 Gaia DR2 references, χ²/ν=0.59, q_plx=h in VizieR). The MJD scale, RV, and correlations are unknown: the CSV leaves them empty, `resolve` returns `ref_epoch_scale=null`, and the model separately states `reference_epoch_scale_assumed_utc`, `rv_assumed_zero`, and the unpublished RA-error basis. These are computational assumptions, not catalogue additions.

The preset propagates this position, uses a geometric geocentric prediction, and centres at `center_mjd=61358`. The number 61358 is the field-centre date, **not WISE's reference epoch or each exposure's date**. The exact specified pair **RA 08h55m10.831s, Dec −7d14m42.53s** converts to **133.79512916666664°, −7.245147222222222°** and remains a fixed centre. No new reference epoch is assigned. In the saved response, `target_info`, `motion`, `center_mjd`, and all `predicted_px` are null.

For coverage, the WISE model was projected **separately**, without changing that response, into its WCS at each date. The origin is the centre of the zero pixel; margin is measured to the pixel-centre bounds `[0,58]` of the 59×59 grid. This is a conservative geometric check, not a check for absent aperture masks.

| UTC date / obs_id | Model RA, Dec, degrees | `[x,y]` in the fixed field FITS grid | Common grid / original cutout margin, px |
| --- | --- | --- | --- |
| 2025-05-05 / `2025W19_1B_0116_1` | 133.761138595, −7.242275669 | [48.7382, 30.6802] | 9.2618 / 24.2235 |
| 2025-11-10 / `2025W46_1A_0075_4` | 133.760199380, −7.242248785 | [49.2836, 30.6959] | 8.7164 / 23.5321 |
| 2025-12-12 / `2025W50_2A_0232_1` | 133.759974557, −7.242251472 | [49.4141, 30.6943] | 8.5859 / 24.2755 |
| 2026-04-21 / `2026W17_1B_0227_2` | 133.758958160, −7.242109393 | [50.0044, 30.7774] | 7.9956 / 22.7751 |

All four predictions lie within FOV=6′, including a 7 px geometric margin for the background annulus. The minimum margin in the preset field is 27.1207 px. This result applies to these dates and cutouts; 6′ is not claimed as a universal minimum.

Detection is **not accepted**. In the 2025-11-10 preset the only `ok/fair` candidate is about 104″ from the prediction. On two later dates candidates lie approximately 0.8″ and 1.4″ from it, but responses remain `ambiguous`: independent registration, PSF photometry with errors, neighbour/artifact checks, and significance accounting for masks are missing. Locally clean pixels, SOURCE, and a high-S/N algorithm-selected spot do not establish a faint target. A detection claim requires these checks, comparison across dates with expected displacement, and flux upper limits where needed. Lack of confirmation is also **not a measured non-detection**. WISE/Spitzer detections from [Luhman (2014)](https://arxiv.org/abs/1404.6501) do not automatically transfer to SPHEREx.

### RW Aur and bright-source masks

New SIMBAD entries associate 156430817820015232 with B and 156430822114424576 with A. The literature identification is reversed: [Dodin et al. (2019), note 1](https://arxiv.org/html/1811.00449) call **Gaia DR2** 156430822114424576 component B; [Kurtovic et al. (2024), §3.1 and §4.6](https://arxiv.org/html/2407.18828v1#S3.SS1) place B southwest of A and note RUWE around 16 for A. The Gaia DR3 coordinates put the second source relative to the first at east≈−1.438″, north≈−0.402″, separation≈1.493″. This supports choosing eastern literature A; matching identifiers across releases alone does not replace the geometric check.

At 6.15″/px the pair occupies about 0.24 px and is unresolved here. Centroid and flux depend on A, B, surroundings, band and masks. Core saturation/clipping and different component spectra can change the apparent centre without measuring A's motion. TRANSIENT/SUR_ERROR/BLOOM occur near the prediction in all three frames; OVERFLOW is absent in the first and present in the second and third. This is more precise than saying all three frames are saturated.

All 43 public-final responses have `quality.verified=false`; even `ok/fair` describes the algorithm, not a verified identification. For all eight targets `motion.pm_measured_arcsec_yr=null`; the WISE coordinate query has no motion model. The algorithm's `saturated_core` flag is built from the overall `contam` mask, which includes OVERFLOW/NONLINEAR as well as ghosts, crosstalk and other contamination. Physical saturation is therefore assessed from original FLAGS, not this label alone. SOURCE denotes an expected known source, not a defective pixel or new detection. It is used in background estimation but should not automatically be removed from target flux. [IRSA definitions](https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/sp.html).

### R Dra: dates and band suitability

All 11 dates are distinct, `n_exp=1`, span_days=0; they are not merged into one long exposure. The baseline is **413.9768 days**, with a 203.0812-day gap between 2025-06-06 and 2025-12-26. The table below shows nearby local bands around 1.65 and 2.40 μm, but all 11 neighbourhoods have OVERFLOW/SUR_ERROR/BLOOM and all centroids are `poor`. Even the 1.6553/1.6502 or 2.4032/2.4013 μm pairs are not usable photometry without checking masks, response and PSF. The whole D3 series cannot be treated as a single band. A calibrated time series, sampling-window and alias analysis, and period uncertainty are absent. The preset label “period ~245 days” is reference context; **the period was not measured from this sample**.

### Proxima Cen, 61 Cyg A and Vega

`%I.0`, `%J`, `%P`, `%X`, and bibliographic codes from saved [Proxima](../../data/validation/responses/proxima-cen-sesame.txt), [61 Cyg A](../../data/validation/responses/61-cyg-a-sesame.txt), and [Vega](../../data/validation/responses/vega-sesame.txt) responses were compared with `resolve`. Proxima/61 Cyg A position and PM cite 2020yCat.1350....0G; Vega cites 2007A&A...474..653V. These references do not verify the reference epoch of the complete solution. Sesame documents `%J` as a J2000 position and `%J.E` as errors and source; this does not justify automatically assigning J2000.0 or J2016.0 to all parameters. [CDS Sesame](https://cdsarc.cds.unistra.fr/doc/sesame.htx).

The application correctly retains `ref_epoch/ref_epoch_format/ref_epoch_scale=null`, `coordinate_model=fixed_catalogue`, `reference_epoch_unverified`, and `proper_motion_not_applied`. Retrieved PM values are not declared zero or used for fictitious propagation. Equivalent coordinates:

| Name | Canonical source_id from response | Fixed RA, Dec, degrees | Comparison result |
| --- | --- | --- | --- |
| Proxima Cen | `NAME Proxima Centauri` | 217.42894222, −62.67949019 | same 4 obs_ids and MJDs |
| 61 Cyg A | `*  61 Cyg A` (two spaces retained) | 316.72474829, +38.74941732 | same 5 obs_ids and MJDs |
| Vega | `* alf Lyr` | 279.23473479, +38.78368896 | same 4 obs_ids and MJDs |

Comparison used `view.frames` and `coordinate_view.frames` themselves, not only the CSV `same_observations` field. It proves equivalence of **field queries**, not correctness of a high-proper-motion star's current position. Astrometry requires a consistent solution with a verified epoch and independent source identification. The 5×5 masks for these three targets were examined at the fixed centre; a clean centre region does not establish the quality of the star after it has moved.

## public-final observation log

All rows: FOV=6′, `single`, one exposure. D3 for all except WISE (D6). Complete product URLs, exact MJD-AVG times and counts for each bit are saved in [evidence.json](../../data/validation/science-acceptance/evidence.json). Here MJD is rounded as in the API; λ is the local band in the response, not the detector's full range. For the WISE objective query, local λ refers to the fixed centre and differs slightly from the preset: **4.9684, 4.4163, 4.6315, 4.4257 μm** for the same four obs_ids.

5×5 region bits: T=TRANSIENT, O=OVERFLOW, E=SUR_ERROR, N=NONFUNC, C=COLD, P=PERSIST, B=BLOOM; S is the number of SOURCE pixels out of 25. For the five presets the region is centred on the prediction; for the three Sesame targets, on the fixed catalogue position. This is local diagnostics, not a map of the whole object.

| Target | UTC date / MJD | obs_id | Release | λ, μm | FLAGS 5×5 | Algorithm / quality |
| --- | --- | --- | --- | --- | --- | --- |
| Barnard | 2025-08-15 / 60902.8654 | `2025W33_2C_0321_1` | DR1 | 2.3074 | T,O,E,B; S=25 | ok/poor |
| Barnard | 2026-03-19 / 61118.3099 | `2026W12_1A_0598_1` | DR1 | 1.7935 | T,O,E,B; S=25 | ok/poor |
| Barnard | 2026-04-10 / 61140.0307 | `2026W15_2A_0137_1` | DR1 | 1.7932 | T,O,E,B; S=25 | ok/poor |
| Barnard | 2026-05-04 / 61164.1338 | `2026W19_1A_0046_1` | DR1 | 2.1062 | T,O,E,B; S=25 | ok/poor |
| Barnard | 2026-08-16 / 61268.6315 | `2026W33_2A_0594_4` | QR3 | 1.6504 | T,O,E,B; S=25 | ok/poor |
| Eps Indi B | 2025-04-29 / 60794.0655 | `2025W18_1B_0081_2` | DR1 | 1.6470 | N; S=25 | ambiguous/— |
| Eps Indi B | 2025-05-20 / 60815.1319 | `2025W21_1B_0227_1` | DR1 | 1.7637 | S=25 | ambiguous/— |
| Eps Indi B | 2025-06-10 / 60836.1875 | `2025W24_1A_0245_2` | DR1 | 1.8475 | S=25 | ambiguous/— |
| Eps Indi B | 2025-09-28 / 60946.9648 | `2025W39_2A_0693_3` | DR1 | 1.6383 | S=25 | ambiguous/— |
| Eps Indi B | 2025-11-01 / 60980.4452 | `2025W44_2B_0494_1` | DR1 | 1.7401 | S=25 | ambiguous/— |
| Eps Indi B | 2026-04-29 / 61159.1631 | `2026W18_1A_0444_1` | DR1 | 1.8050 | S=25 | ambiguous/— |
| Eps Indi B | 2026-09-22 / 61305.2276 | `2026W39_1A_0276_1` | QR3 | 1.7215 | S=25 | ambiguous/— |
| WISE 0855 | 2025-05-05 / 60800.7007 | `2025W19_1B_0116_1` | DR1 | 4.9674 | S=15 | ambiguous/— |
| WISE 0855 | 2025-11-10 / 60989.2975 | `2025W46_1A_0075_4` | DR1 | 4.4178 | T,E,N,C; S=14 | ok/fair |
| WISE 0855 | 2025-12-12 / 61021.2835 | `2025W50_2A_0232_1` | DR1 | 4.6336 | S=19 | ambiguous/— |
| WISE 0855 | 2026-04-21 / 61151.1986 | `2026W17_1B_0227_2` | DR1 | 4.4237 | S=18 | ambiguous/— |
| RW Aur | 2025-09-11 / 60929.6781 | `2025W37_2A_0051_1` | DR1 | 2.3540 | T,E,B; S=25 | ok/poor |
| RW Aur | 2026-02-20 / 61091.7203 | `2026W08_2A_0233_1` | DR1 | 1.9367 | T,O,E,B; S=25 | ok/fair |
| RW Aur | 2026-09-13 / 61296.2983 | `2026W37_2A_0581_1` | QR3 | 2.1782 | T,O,P,B; S=25 | ok/fair |
| R Dra | 2025-04-24 / 60789.8803 | `2025W17_4B_0021_1` | DR1 | 1.6627 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2025-05-16 / 60811.9565 | `2025W20_2D_0243_1` | DR1 | 1.6553 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2025-06-06 / 60832.1375 | `2025W23_2B_0129_2` | DR1 | 2.4032 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2025-12-26 / 61035.2187 | `2025W52_2A_0161_2` | DR1 | 1.6419 | O,E,B; S=25 | ok/poor |
| R Dra | 2026-02-08 / 61079.4825 | `2026W06_2A_0609_1` | DR1 | 2.4013 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2026-03-01 / 61100.4643 | `2026W09_2A_0550_1` | DR1 | 1.8241 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2026-03-21 / 61120.6976 | `2026W12_2A_0398_2` | DR1 | 1.6502 | T,O,E,P,B; S=25 | ok/poor |
| R Dra | 2026-04-11 / 61141.4728 | `2026W15_2A_0483_1` | DR1 | 1.9202 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2026-05-02 / 61162.0467 | `2026W18_2A_0354_3` | DR1 | 1.6297 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2026-05-22 / 61182.1332 | `2026W21_2A_0104_1` | QR3 | 1.7610 | T,O,E,B; S=25 | ok/poor |
| R Dra | 2026-06-12 / 61203.8571 | `2026W24_2B_0242_1` | QR3 | 1.8251 | T,O,E,B; S=25 | ok/poor |
| Proxima | 2025-07-22 / 60878.4934 | `2025W30_1B_0271_1` | DR1 | 1.7134 | S=25 | ambiguous/— |
| Proxima | 2025-08-23 / 60910.2220 | `2025W34_2B_0324_1` | DR1 | 1.7132 | S=25 | ambiguous/— |
| Proxima | 2026-02-15 / 61086.9318 | `2026W07_2A_0639_2` | DR1 | 1.6564 | P; S=25 | ambiguous/— |
| Proxima | 2026-07-20 / 61241.3548 | `2026W30_1B_0070_1` | QR3 | 1.6960 | S=25 | ambiguous/— |
| 61 Cyg A | 2025-05-25 / 60820.1708 | `2025W21_2C_0442_1` | DR1 | 1.7070 | S=25 | ambiguous/— |
| 61 Cyg A | 2025-06-16 / 60842.0471 | `2025W25_1B_0010_1` | DR1 | 2.2257 | S=25 | ambiguous/— |
| 61 Cyg A | 2025-10-24 / 60972.8218 | `2025W43_2A_0295_1` | DR1 | 1.8175 | N; S=25 | ambiguous/— |
| 61 Cyg A | 2026-05-25 / 61185.6566 | `2026W22_1A_0106_1` | QR3 | 1.9757 | S=25 | ambiguous/— |
| 61 Cyg A | 2026-07-05 / 61226.5961 | `2026W27_2A_0574_1` | QR3 | 1.6831 | S=25 | ambiguous/— |
| Vega | 2025-08-28 / 60915.5635 | `2025W35_2A_0014_1` | DR1 | 1.9103 | O,E,N,B; S=25 | ok/poor |
| Vega | 2026-04-05 / 61135.2197 | `2026W14_2A_0517_1` | DR1 | 2.1497 | O,E,B; S=25 | ok/poor |
| Vega | 2026-05-04 / 61164.4098 | `2026W19_1A_0141_1` | DR1 | 2.1493 | O,E,P,B; S=25 | ok/poor |
| Vega | 2026-08-30 / 61282.4707 | `2026W35_2A_0569_1` | QR3 | 1.8796 | O,E,B; S=25 | ok/poor |

The exact WISE query has the same four obs_ids and dates; all four algorithm results are `ambiguous`. `quality` and `wavelength_samples` were not saved in that older response: these are unknown values, not `verified=false` or zero errors carried over from the newer preset. An additional check of its real FITS found the same local bits at the model position (including T/E/N/C on 2025-11-10).

## Geometry and model limits for D3

All 43 entries in [geometry_audit.json](../../data/validation/geometry_audit.json) matched public-final by target, obs_id and date. The original audit found minimum common-grid margin **26.9559 px**, original-cutout margin **36.2833 px**, and maximum PNG marker discrepancy from WCS recalculation **0.000488 px**. Recalculation from corresponding available FITS confirmed coverage and the correct single y flip. This checks calculations with the same predictor, **not independent registration on stars**. The exact-coordinate WISE query was not in those 43 entries; its four cutouts were checked separately above, and their smaller margin must not be replaced by the overall minimum.

The five-target model propagates catalogue position in space and applies geometric geocentric parallax in ICRS axes. It does not include aberration, light deflection, spacecraft position in place of the geocentre, or a new covariance at each date; Ba/Bb orbital propagation is also absent. Formal catalogue errors do not assess masks/PSF/blending. Aligning an absolute marker with a visible source requires empirical registration and checking the particular WCS conventions.

## Specific changes and reasons

1. **OBJECT_PARAMETERS.csv was retained unchanged.** All five rows agree with their sources; missing RV, time scales and correlations were not replaced with zeros. The original 15 predictions and initial Barnard series were not overwritten.
2. In the Epsilon Indi B row of **OBJECT_PARAMETERS.md**, the HTML Table 2 link was replaced by Table 2 in the PDF, p. 3790: the previous HTML URL lacks the required lower part of the table. The Ba identification itself is correct; the change provides access to the precise primary evidence.
3. In the WISE row of **api/resolve.py / PRESETS**, “visible only at 4–5 µm” was changed to “D6: 4.42–5.00 µm; detection requires verification.” Detector choice does not prove visibility in SPHEREx, and “only” unjustifiably restricted the object's physical observability. Evidence: D6 range, original [WISE/Spitzer observations](https://arxiv.org/abs/1404.6501), and ambiguous application responses. Field settings and computational model did not change.
4. The RW Aur row in **TARGET_VALIDATION.md** was clarified against actual FLAGS: masks are present on every date, OVERFLOW on two. Links to this review and its limits were added to TARGET_VALIDATION.md and DILAYLA_ACCEPTANCE.md.

## Acceptance of D2, D3 and D6

| Item | Accepted by this review | Requires additional data |
| --- | --- | --- |
| **D2.1–4** | Sources of five rows, Gaia IDs as strings, Ba and eastern literature RW Aur A components; MJD 57633.39 for WISE; canonical IDs of three Sesame stars; original responses and provenance retained | A more precise Ba/Bb system solution needs an absolute reference position with epoch and orbit; RW Aur's letter ambiguity must remain visible |
| **D2.5** | Unknown parameters remain empty/null; assumed RV=0 and UTC in the model are explicitly identified | Verified WISE MJD scale and covariance; complete-solution epochs for the three Sesame stars if propagation is required |
| **D2.2, D2.6** | Canonical names and real responses for the three stars accepted; encoding/service-error evidence is accepted from the existing technical review | This scientific review is not a new live test of resolver availability |
| **D3.1–4** | Consistent model for five targets; signs of the Barnard prediction; fixed raw RA/Dec; unknown Sesame epoch not invented; actual MJD-AVG UTC in checked FITS | Absolute apparent astrometry needs direction-frame/WCS reconciliation and empirical registration; the three extra stars require verified solutions with epochs |
| **D3.5–6** | Geometric coverage of 43 saved frames, WCS/PNG flip and margin; separately four WISE objective cutouts with ≥7 px margin | Aperture/annulus suitability by masks and background, registration/PSF errors; geometric margin alone does not establish these |
| **D6: five NASA + three additional** | Selected observations identified; distinct dates confirmed; API/browser viewing accepted from saved evidence; names and coordinates gave the same obs_ids and times for the three additional stars | Claims of each target's detection, measured motion or variability require reliable association, registration, comparable bands and measurements with errors. These claims are not accepted |
| **D6: WISE objective** | Exact coordinate pair checked separately from the preset; four D6/DR1 products available in the saved set; prediction inside all 6′ fields | WISE detection and its statistical significance are not confirmed; the checks above are required |
| **D6: R Dra** | Eleven distinct dates and band/mask limitations confirmed | Period and amplitude were not measured from this series; usable calibrated photometry in matched bands and adequate temporal sampling are required |

The final D6 acceptance scope is **reproducible viewing of real fields and an appropriately qualified catalogue prediction**. Full scientific certification B1, agreement by shared-file owners, and permanent hosting are not closed by this report. Timeouts of new IRSA requests are not interpreted as absence of data.

## Verification of changes

`uv run --offline --locked pytest -q`: **197 passed, 1 skipped** (a separate live test was skipped; no new live application pass is claimed). Additional checks covered JSON/hash integrity, all 43+4 entries against original responses, three name/coordinate query comparisons, equality of parameters to source catalogues, preservation of OBJECT_PARAMETERS.csv and old predictions, and `git diff --check`. No new detections or photometric measurements were performed.
