# SPHEREx Sky Viewer — pitch text (draft 1)

Slide text is what goes on the slide. *Notes* are what the speaker says.
Structure follows the jury criteria: problem → solution → validity → creativity → value.
Target length: 4–5 minutes + live demo.

---

## 1. Title

**SPHEREx Sky Viewer**
Blink the infrared sky — see what moves and what changes.

NASA Space Apps Challenge 2026 · team name · city

*Notes:* One sentence: we built a public web tool that turns SPHEREx archive images into a short movie of any patch of sky.

---

## 2. The problem

- SPHEREx maps the **whole sky every 6 months** in **102 infrared bands**
- More than **a billion objects**: brown dwarfs, asteroids, comets, variable stars
- Discoveries come from **comparing images taken at different times** — that is how Pluto was found, and how volunteers on Backyard Worlds find brown dwarfs
- But the raw data is not public-friendly: 2040×2040 FITS files, 6 detectors, a wavelength that changes across every image, dozens of quality flags

**Too much data for any one person — and too hard for most people to open.**

*Notes:* The data is open, but only experts can use it. A student or a volunteer cannot download 70 MB files, subtract the zodiacal light and line up images by hand.

---

## 3. Our solution

Type a **name** or **RA / Dec** → get a **blinking movie** of that patch of sky from every SPHEREx visit.

- Live data from **NASA/IPAC IRSA** (DR1 + QR3) — nothing is pre-downloaded, so new releases appear automatically
- Choose **detector (1–6)** and **field of view (1–20′)**
- Every view has a **shareable link**: coordinates live in the URL, so the tool can be embedded anywhere
- Works in any browser, no installation

*Notes:* The first request for a new field takes a few minutes because we fetch cutouts from IRSA; after that the result is cached and opens instantly.

---

## 4. How it works

1. **Find** every SPHEREx exposure covering the point (IRSA TAP, with SIA2 as fallback)
2. **Cut out** only the field you asked for — not the full 70 MB file
3. **Clean**: subtract the zodiacal-light model, mask bad pixels and cosmic rays *before* interpolation
4. **Group** exposures into epochs and **co-add** them on one north-up grid
5. **Render** frames, a contamination-flag layer and a JSON with every measurement

*Notes:* Masking before interpolation matters: one hot pixel spread over 2×2 grid pixels looked like a star with signal-to-noise 279. We caught it and fixed the order.

---

## 5. Four ways to see change

| Mode | What you see |
|---|---|
| **Single** | One raw exposure per epoch |
| **Co-add** | All exposures of an epoch stacked — deeper image, less noise |
| **Difference** | Each epoch minus the median of the *other* epochs: a moving object leaves a **+ / −** pair |
| **Trail** | All epochs in **one picture**, coloured by time: blue = first, yellow = last |

*Notes:* Trail is the fastest way to spot motion: a moving star becomes a blue-to-yellow streak, while stars that do not move disappear.

---

## 6. More than pictures

- **Spectrum on the fly**: SPHEREx's filter sees each object at a different wavelength in every exposure — we collect those points into a spectrum across all 6 detectors
- **Light curve** in a fixed wavelength window, so brightness changes are not confused with the shape of the spectrum
- **Contamination flags** overlay: saturation, ghosts, crosstalk, streaks
- **Sky map** from Aladin Lite (CDS) next to the player
- The same engine as a **Python library** (`spxhack`): researchers get numpy arrays and FITS cubes for DS9

*Notes:* The spectrum is built from the same exposures as the movie — no extra observations needed.

---

## 7. Tested on the five challenge targets

| Target | SPHEREx epochs | What the tool shows |
|---|---|---|
| **Barnard's Star** | 5 dates, Aug 2025 – Aug 2026 | ~10″/yr motion, almost due north, about 1.7 pixels per year |
| **Epsilon Indi B** | 7 dates, Apr 2025 – Sep 2026 | The brown-dwarf pair field; the pair is unresolved at 6″ pixels |
| **WISE 0855−0714** | 4 dates, D6 (4.4–5.0 μm) | Field and predicted track of the coldest known brown dwarf |
| **RW Aurigae** | 3 dates, Sep 2025 – Sep 2026 | Young variable star; A and B components are 1.5″ apart (unresolved) |
| **R Draconis** | 11 dates, Apr 2025 – Jun 2026 | Mira variable: a real time series for brightness changes |

Every identity checked against Gaia DR3, SIMBAD and the literature; every observation listed by ID.

*Notes:* We tested the tool on the challenge targets and documented each one: which catalogue record, which exposures, and what can and cannot be concluded.

---

## 8. Science we can stand behind

- Every number comes with its **method, error and quality flag**
- When the data cannot answer, the tool **says so instead of guessing**:
  - saturated core → no fake centre
  - two candidates → both shown, none chosen
  - different wavelengths between epochs → warning, not "variability"
- **About 290 automated tests**, plus replay of saved archive answers
- All checks are documented in the repository

*Notes:* Example: Barnard's Star is so bright that its core saturates. Many tools would still return a position. We show the predicted position and mark the measurement as unreliable. For a tool used by the public, honest limits are part of the science.

---

## 9. What is new in our approach

- **Time-coloured trail**: motion in one still image — no need to watch a video
- **Spectrum assembled from the moving filter**: we turn SPHEREx's changing wavelength into a feature
- **Leave-one-out difference**: the epoch is never compared with itself, so moving objects give a clean + / − pair
- **One engine, two interfaces**: a volunteer spots something on the website, a scientist reproduces it in Python in three lines

*Notes:* The key idea: what a volunteer sees on the website and what a scientist computes in a notebook come from the same code.

---

## 10. Who it is for

- **Volunteers** in projects like Backyard Worlds: a new data set to search for brown dwarfs and moving objects
- **Students and teachers**: watch a star move and a Mira change brightness, using real NASA data
- **Researchers**: a quick look before downloading full files, and a direct path to numpy / FITS
- **Other websites**: any view can be embedded by URL

*Notes:* Most known brown dwarfs and comets were found by volunteers. SPHEREx adds infrared colours that are especially sensitive to cold objects.

---

## 11. What comes next

- Gaia DR3 and WISE catalogue overlays on the sky map
- New SPHEREx releases picked up from configuration, without code changes
- PSF fitting and deblending for crowded fields
- "Flag a candidate" button so volunteers can report what they find

*Notes:* The pipeline is already live against the archive, so each new release makes the tool better without extra work.

---

## 12. Thank you

**SPHEREx Sky Viewer**
live demo link · GitHub: github.com/dmitryceohris/nasa_space_apps2026

Team: names and roles

Data: NASA SPHEREx via NASA/IPAC IRSA · Gaia DR3 (ESA) · SIMBAD and Aladin Lite (CDS)

---

## Live demo script (about 90 seconds)

1. Open a **shareable link** to Barnard's Star in **Trail** mode — show the blue-to-yellow streak.
2. Switch to **Co-add** and press **Play** — the star shifts between epochs; point at the predicted-position marker.
3. Turn on the **flags layer** — show the saturated core and explain why the tool does not invent a centre.
4. Open the **Spectrum** tab — the points come from the same exposures.
5. Switch to **R Draconis** — 11 dates; open the **light curve** in one wavelength window.

Before the demo: run the server locally with a pre-warmed cache, so the demo does not depend on the IRSA connection.
