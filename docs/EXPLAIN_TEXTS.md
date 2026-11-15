# "Explain this region" — texts (draft 1, for approval)

One fixed text per block: the same region always gives the same explanation.
Variety comes from the data (numbers, caveats), not from chance.

Notation:
- `{name}` — placeholder, filled from the measurement (formats below)
- `[term]` — underlined term with a tooltip from the glossary (first occurrence in a card only)

---

## 1. Card layout

```
<Verdict title>
<What we see>
<Numbers>
<Interpretation>            (only if the verdict has one)
<Caveats>                   (0–3, by conditions, in the order listed)
<Next step>
<Footer>                    (always)
```

If the area holds more than one source, the cards are preceded by:

> We found {n_sources} sources in this area. Each one is explained separately below.

---

## 2. Verdicts

Checked in this order; the first that fits is shown.

Implementation note: a **saturated** source is checked before "Likely artefact". On real data Barnard's Star sits
on BLOOM pixels (charge spilling from its saturated core) and would otherwise be called "not a real object".

### 2.1 Not enough data

Condition: fewer than 2 [epochs] with data in the area.

> **Not enough data**
> This area has usable data in only {n_epochs_with_data} [epoch]. At least two are needed to see any change.

Next step: *Try a larger field of view or another detector.*

### 2.2 Likely artefact

Condition: the source sits on pixels with contamination flags in the epoch where it appears.

> **Likely artefact**
> This is probably not a real object. The pixels here carry the {flag_names} [flags], which SPHEREx uses to mark instrument effects rather than light from the sky.

Next step: *Turn on the flags layer to see the affected pixels.*

### 2.3 Can't tell

Condition: saturated core, or two candidates and no way to choose, or the source touches the edge of the data.

> **Can't tell**
> {reason_text}

`{reason_text}` — one of:
- saturated: The centre of this source is [saturated]: it is so bright that the detector reached its limit, so its exact position and brightness cannot be measured.
- ambiguous: Two sources are too close to tell which one changes.
- edge: The source touches the edge of the data, so its centre would be biased.

Next step: *Look at the Co-add view: the source is real, only the measurement is limited.*

### 2.4 Moving object candidate

Condition: a [residual] pair (bright + dark spot) in at least 2 epochs, and the centre shifts consistently between epochs.

> **Moving object candidate**
> Something here is on the move. In {n_pair} of {n_epochs} epochs a bright spot appears in a new place and leaves a dark gap where it used to be.
> Between {date_first} and {date_last} it shifted by {shift_arcsec}″ ({shift_px} [pixels]), about {pm_arcsec_yr} ± {pm_err}″ per year toward the {direction}.

Interpretation — one of:
- matches prediction: This agrees with the predicted path of {target} to within {offset_px} px.
- no prediction, fast (≥ 0.5″/yr): Objects that drift this fast are usually nearby stars or [brown dwarfs]: the closer an object is, the faster it seems to move. This is a measured [proper motion] candidate, not a confirmed one.
- no prediction, slow (< 0.5″/yr): The shift is small. It may be real motion, or the images may be lined up slightly differently.
- disagrees with prediction: The shift does not match the predicted path of {target} (off by {offset_px} px). This may be a different object.

Next step: *Open the Trail view: the motion shows as a blue-to-yellow streak.*

### 2.5 Seen only once

Condition: a bright residual in exactly 1 epoch and nothing at that place in the others.

> **Seen only once**
> A source appears here on {date_seen} and is absent in the other {n_other} epochs.
> Fast movers such as [asteroids] and comets cross this field in days, so they show up in just one epoch. A [cosmic ray] or an instrument effect can look the same.

Next step: *Check this date in the Single view: a real object has the shape of a star, a cosmic ray is one sharp pixel.*

### 2.6 Brightness change

Condition: a residual of one sign at a fixed place; no shift.

> **Brightness change**
> The source stays in place, but its brightness changes by about {change_pct}%: brightest on {date_bright}, faintest on {date_faint}.
> This is how [variable stars] look in the data.

Next step: *Open the Light curve tab to follow the brightness at one wavelength.*

### 2.7 Nothing significant

Condition: no residual above 5 [S/N] in the area.

> **Nothing significant**
> Nothing in this area changes beyond the [noise] level across {n_epochs} epochs.

Next step: *Select a smaller area around a source, or try the Difference view.*

---

## 3. Caveats

Added after the interpretation, in this order, when the condition holds.

| Condition | Text |
|---|---|
| wavelength range between epochs is wider than one channel | The epochs were observed at different wavelengths ({wave_min}–{wave_max} µm). Part of the change may come from the object's [spectrum], not from real [variability]. |
| moving, shift < 2 px | The shift is smaller than two pixels, so the bright and dark spots partly overlap: treat the speed as approximate. |
| moving, known [parallax] gives ≥ 20% of the shift | Part of this shift is [parallax]: a yearly wobble caused by the Earth moving around the Sun. |
| moving or brightness, some epochs have saturated or flagged pixels on the source | In {n_flagged} epochs the source has flagged pixels; those epochs were not used for the numbers above. |
| data coverage in the area < 90% in some epoch | Part of this area has no data in {n_partial} epochs. |
| registration not verified (always for moving) | Image alignment was taken from the telescope's coordinates and was not checked against reference stars. |

---

## 4. Footer (always)

> This is an automatic reading of the data, not a confirmed discovery.

---

## 5. Placeholder formats

| Placeholder | Format | Example |
|---|---|---|
| `{n_sources}`, `{n_epochs}`, `{n_pair}`, `{n_other}`, `{n_flagged}`, `{n_partial}`, `{n_epochs_with_data}` | integer | 5 |
| `{date_first}`, `{date_last}`, `{date_seen}`, `{date_bright}`, `{date_faint}` | `MMM D, YYYY` | Aug 15, 2025 |
| `{shift_arcsec}`, `{shift_px}`, `{offset_px}` | 1 decimal | 10.6 |
| `{pm_arcsec_yr}`, `{pm_err}` | 1 decimal | 10.3 ± 0.8 |
| `{direction}` | one of 8: north, north-east, east, south-east, south, south-west, west, north-west | north |
| `{change_pct}` | integer, rounded to 5 | 35 |
| `{wave_min}`, `{wave_max}` | 2 decimals | 1.63–2.40 |
| `{target}` | display name | Barnard's Star |
| `{flag_names}` | human names joined with "and" (see 6) | ghost and crosstalk |

Plural forms: "1 epoch" / "2 epochs", "1 pixel" / "2 pixels".

---

## 6. Flag names shown to the user

| Flag | Shown as |
|---|---|
| GHOST, GHOST_FPA, GHOST_EXT | ghost |
| CROSSTALK | crosstalk |
| STREAK | streak |
| BLOOM | bloom |
| SNOWBALL | snowball |
| HALO | halo |
| OVERFLOW | saturation |
| NONLINEAR | nonlinear response |

---

## 7. Glossary (tooltips)

| Term | Tooltip |
|---|---|
| epoch | A group of SPHEREx images of the same place taken within a few days. SPHEREx returns to each part of the sky about every six months. |
| pixel | The smallest element of the image. One SPHEREx pixel covers 6.15″ of sky. |
| arcsecond (″) | A unit of angle: 1/3600 of a degree. The full Moon is about 1800″ across. |
| proper motion | How fast a star drifts across the sky, in arcseconds per year. Nearby stars drift fastest. |
| parallax | A small yearly back-and-forth shift of nearby objects, caused by the Earth's orbit around the Sun. |
| residual | What is left after subtracting the other epochs. Bright means more light now, dark means less. |
| saturated | So bright that the detector reached its limit, so the true brightness and centre are unknown. |
| flags | Labels that SPHEREx attaches to pixels affected by instrument problems. |
| ghost | A faint false image of a bright star, caused by reflections inside the telescope. |
| crosstalk | A false signal copied from a bright source into another part of the detector by the electronics. |
| cosmic ray | A charged particle that hits the detector and lights up one or a few pixels for a single image. |
| noise | Random variation in every image. Changes smaller than the noise cannot be told apart from chance. |
| S/N | Signal-to-noise ratio: how many times a signal is stronger than the noise. Below 5 we do not trust it. |
| spectrum | How bright an object is at each wavelength of light. |
| wavelength | The "colour" of light. SPHEREx sees infrared light from 0.75 to 5.0 µm, invisible to the eye. |
| variability | A real change in the object's own brightness over time. |
| variable stars | Stars whose brightness changes, for example pulsating Mira stars or young stars still gathering gas. |
| brown dwarfs | Objects heavier than planets but too light to shine like stars. They are cold and glow mainly in infrared light. |
| asteroids | Small rocky bodies orbiting the Sun. Seen from Earth they move very fast across the sky. |
| detector | SPHEREx has six detectors; each covers its own range of infrared wavelengths. |
