# Scenario validation

## Combined acceptance, November 14, 2026

Current checked code **be2ca1a**, base **afa35ef**, scientific port **cb29e02**. [processing-acceptance.md](reviews/processing-acceptance.md) concludes B1 with **“accepted with limitations”**: field viewing and catalogue predictions are accepted; independent registration and measurement accuracy are not confirmed. A new replay of saved real FITS used 47 inputs, 36 /api/view scenarios, 188 frames, and 0 archive data requests; [responses](../data/validation/processing-acceptance/replay.json). Full suite: **259 passed, 1 skipped**; `git diff --check` clean.

The information below records live and browser checks at the specified commits. The FITS replay is not a new live external run; controlled browser checks were separately repeated in processing-acceptance/browser.json. D2–D7 are closed for viewing according to the [acceptance map](DILAYLA_ACCEPTANCE.md); permanent hosting B3 remains separate. Direct participant confirmation of processing/observations was not obtained.

November 14, 2026. Checked code: `aa7ce46`, base `39bdbd5`. Current report changes are identified separately in Git. D1 was already in the base. D2–D5 were refined in this run; D6 was checked to the extent described below. Shared changes were already merged into the base; no new agreement with participants is claimed.

## Real data

The server ran with `PREWARM=0`, without mock data. The check uses `single` mode, a 6′ field, one exposure per epoch, D6 for WISE and D3 for the other objects. Manifest: [target_cases.csv](../data/validation/target_cases.csv). Responses: [responses](../data/validation/responses/). The JSON records source, components, actual reference epochs, assumptions, obs_ids, times and warnings. Raw Sesame responses for three additional stars are stored separately (`*-sesame.txt`).

| Target | Frames on different dates | Limitation |
| --- | ---: | --- |
| Barnard | 5 | Contamination/saturation, poor centroid quality, different bands |
| Epsilon Indi B | 7 | Ba coordinates; ambiguous candidates, Ba/Bb unresolved |
| WISE 0855 | 4 | Ambiguous candidates, different bands; target detection not confirmed |
| RW Aur | 3 | Eastern literature A; A/B unresolved; TRANSIENT/SUR_ERROR/BLOOM in all three prediction neighbourhoods, OVERFLOW in two; see scientific review below |
| R Dra | 11 | Saturation and different bands; period not measured |
| Proxima Cen | 4 | Sesame J2000 without a verified reference epoch; fixed position, ambiguous candidates |
| 61 Cyg A | 5 | Same epoch limitations; different bands and ambiguous candidates |
| Vega | 4 | Same epoch limitations; saturation and poor centroids |

The exact WISE pair `08h55m10.831s, -7d14m42.53s` separately returned 4 frames at FOV=6′, D6; [response](../data/validation/responses/wise-objective.json). The centre remained fixed, `target_info=null`, and no motion prediction was added.

The first pass made real HTTP requests and FITS downloads from IRSA; some metadata and FITS may already have been in the local cache. R Dra exceeded the initial 90-second limit; the job completed later and returned 11 frames. All eight targets were rerun on the final version. It used real FITS and cache; this is **not** evidence of a cold start without cache or current completeness of the archive sample. IRSA failures/503 are not interpreted as absence of coverage.

## Browser

The main HTTP page was checked for all eight objects: PNGs loaded (59×59 pixels), and the navigation button changed the first epoch date to the second. Log: [browser-real.txt](../data/validation/browser-real.txt). Saved links with `target`, `fov`, `det`, `mode`, and `max_per_epoch` were checked. At 375 px width the image occupies 351 px; warnings persist after advancing to the next frame. After the log was recorded only the centroid label was clarified to “centroid candidate”.

Controlled browser checks:

- The empty form requests only `/api/targets`.
- `target=barnard&ra=133.795` and pressing Light curve show a prompt; no data requests occur.
- When RA changes quickly from 10 → 30 and the RA=10 response is delayed, the RA=30 warning remains.
- A PNG failure shows an error, the main image has no `src`, and plots do not start automatically.
- Images work when Plotly is unavailable; plots require a separate click.

New server tests cover outgoing encoding of special characters in names, resolver unavailability, provenance/unknown values, five targets offline, agreement with `sky_motion` on three dates, UTC/TAI and exposure-midpoint priority, different missing-data states, and the absence of a reliable marker for a poor centroid.

## Reproduction

From the repository root:

```sh
PREWARM=0 .venv/bin/python -m uvicorn api.main:app --host 127.0.0.1 --port 8011
.venv/bin/python scripts/validate_targets.py --images --timeout 150
.venv/bin/python -m pytest -q
```

Links in the CSV are relative: open them at the running application's address. The script writes a new manifest and responses for the selected run; FITS files are not added to Git. Last full test run: **195 passed, 1 skipped**. The skipped test was the separate live check `test_live_barnard_cutout`; real API/browser requests are described above.

## Continued validation: actual result

A cold pass used an empty separate API/FITS cache, disabled local tables, and `CENTER_MJD=61358`: [cold/target_cases.csv](../data/validation/cold/target_cases.csv). All eight stars returned real frames. Equivalent coordinate queries returned the same obs_ids for Proxima Cen, 61 Cyg A and Vega. A repeat pass, Linux container with persistent volume, and requests through a temporary external address were then checked: [repeat](../data/validation/repeat/target_cases.csv), [container](../data/validation/container/target_cases.csv), [public-preview](../data/validation/public-preview/target_cases.csv).

[geometry_audit.json](../data/validation/geometry_audit.json) contains 43 frames: the axis/single y-flip check passed; the target position has a geometric margin of 7 pixels for the background annulus on the common grid and original FITS. This is not an independent measurement of registration using reference stars.

Browser checks are saved as repeatable scenarios: `scripts/check_viewer_browser.js` covers controlled errors, late responses, name/mode/tab changes, candidates, mobile width, and Plotly/PNG failures; `scripts/check_targets_browser.js` covers real images of eight stars and date switching. Logs: [browser-controlled.txt](../data/validation/browser-controlled.txt), [browser-public.txt](../data/validation/browser-public.txt).

HTTP checks also covered one epoch in all four modes, absence of a false measurement of a stationary star when background/NaN changes near an edge, and recovery of a partial response after a short cache lifetime. Metadata refresh, job completion race, queue limit, and maintaining one job until already-started network calls finish after timeout were checked separately.

The container was built and started, and its healthcheck passed. The persistent volume survived container replacement. Instructions: [DEPLOYMENT.md](DEPLOYMENT.md).

## Open criteria

- Direct team confirmation of selected observations/components and new contract fields has not been obtained.
- B1 is now accepted for viewing with limitations (be2ca1a); empirical registration and scientific certification of measurements remain open for quantitative claims. Centroids are labelled as candidates; quality is not declared verified; the difference scale is not declared calibrated statistical significance.
- B3: no permanent server/domain was selected (the user said “none”). The temporary address below depends on containers running on this computer. Testing from a separate physical device had not yet been confirmed.

Temporary demo: https://three-employers-look-obtaining.trycloudflare.com/ . Its availability does not imply permanent hosting. Quick tunnels are intended for testing: [Cloudflare documentation](https://developers.cloudflare.com/tunnel/get-started/quick-tunnels/). Old scientific tables and predictions were not overwritten.

### Partial result snapshot for push

The latest repeated external pass was stopped at the request to save and push the current state. The first three responses in `public-preview/responses` were already updated; the rest belong to the previous successful pass. The latest browser attempt through the external address exceeded the initial 30-second limit; a retry with a 180-second limit was not completed. `browser-public.txt` records the incomplete check and is not evidence that all eight targets passed browser checks at the external address. The earlier local browser check remains separate evidence.

### Completed external pass after push

A new complete pass on code `aa7ce46` finished: [public-final/target_cases.csv](../data/validation/public-final/target_cases.csv). All eight results have `available`; three additional stars yield the same observations with equivalent coordinates. The browser separately loaded real PNGs for all eight targets and switched the first date to the second. Structured result: [browser-public-final.json](../data/validation/browser-public-final.json), [log](../data/validation/browser-public-final.txt). This pass supersedes the incomplete attempt as external-check evidence. It is an external HTTP address opened from the work computer; a separate physical device and independent permanent server were not checked by this pass.

### Final route check

The external check found a 404 at `/dev` without a slash. An explicit redirect to `/dev/` was added and tested after rebuilding the container. Result: `/dev` → 307, `/dev/` → 200; invalid and incomplete coordinates → 400; PNGs for all eight targets → 200 with a valid PNG signature. Evidence: [public-routes.json](../data/validation/public-routes.json). Last full suite: **196 passed, 1 skipped**. The running container's source matches this version; the volume with real FITS was preserved during replacement.

The D7 handoff is prepared: contract in API.md, startup in DEPLOYMENT.md, observations and browser evidence in data/validation. B3 remains open pending permanent hosting and testing from another device. The user was asked about hosting and phone testing; no response had yet been received.

### Feedback from testing on other devices

The user reported that one tester confirmed it worked, while a second could not find the date. This is user confirmation of operation on a separate device; no separate device log was supplied. In response, the date was moved to a large separate line, date controls placed directly below it, and buttons labelled. Permanent hosting remains unconfirmed.

After updating the external demo, a browser check showed the large Barnard date changing from 2025-08-15 to 2026-03-19, with the counter Date 2 of 5. At 375 px width the date uses a 21 px font; buttons are immediately adjacent. Log: data/validation/browser-date-visibility.txt.

Detailed acceptance map for each item: [acceptance map](DILAYLA_ACCEPTANCE.md). The latest check adds an HTTP scenario for an out-of-frame prediction: **197 passed, 1 skipped**; application code was not changed after the date-visibility change.

### Scientific review of sources and selected observations

[science-acceptance.md](reviews/science-acceptance.md) compares all 43 public-final frames and four frames from the exact WISE coordinate query against original products, epochs, FITS masks, and primary publications. This reconciles identity and interpretation limits for D2, D3 and D6; detection of every target, measured motion, and the R Dra period are not confirmed. The application flag `saturated_core` combines several kinds of contamination and does not by itself prove saturation of a particular target. Release information was recovered from saved URIs and available FITS; new IRSA requests timed out during this review, and current archive completeness was not checked.
