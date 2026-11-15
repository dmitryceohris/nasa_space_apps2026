# SPHEREx Blink: user guide and interpreting results

SPHEREx Blink lets you explore SPHEREx observations at different times. It is an exploratory viewer: a visible spot, a model prediction or a plotted point alone is not proof of a source detection, motion or variability.

## Start with a target

1. Open the [hosted viewer](https://spherex-sky-viewer-zeta.vercel.app/).
2. Select a **Test object**, enter an **Object name**, or provide both **RA** and **Dec**. Coordinates take priority if both coordinates and a name are present. RA and Dec can be decimal degrees; sexagesimal RA is in hours and sexagesimal Dec is in degrees (for example, `08h55m10.831s` and `-7d14m42.53s`).
3. Choose a **Field of view** from 1 to 20 arcminutes and a **Detector** from 1 to 6. The available wavelength coverage and observations depend on the detector and field.
4. Choose **Image source**. **Latest IRSA** searches available archive collections and processes observations; it may take time, and new observations can appear later than in other archive services. **Saved observations** opens a fixed, previously processed selection when one exists for the target or exact supported coordinates. Saved fields have fixed settings; if no saved selection matches, use Latest IRSA.
5. Wait for the images to load. If no frames are returned or an error appears, read the message and try another source or selection. A missing result does not establish that the object is absent from the archive.

On the hosted viewer, changing the field of view or detector switches the image source to Latest IRSA. Moving the slider updates its displayed value immediately; loading starts when the change is completed. Returning to saved observations restores the field and detector of a matching saved selection. The GitHub checkout can differ from the hosted version.

## Compare observations

Choose one of the image modes:

| Mode | What it shows | How to read it |
| --- | --- | --- |
| **single** | One selected exposure per epoch. | Useful for inspecting individual observations; coverage and noise can vary. |
| **coadd** | Exposures grouped and combined by epoch on a common grid. | A frame can contain multiple exposures, possibly over a range of dates and wavelengths. |
| **diff** | Each epoch compared with a reference made from the other epochs. | Coloured positive and negative residuals show relative differences. The residual scale is not calibrated statistical significance. |
| **trail** | A time-coloured visualisation of residuals across epochs. | Colour indicates time, not a verified trajectory. This is a combined image, not a series of independent position measurements. |

Difference and trail views need multiple epochs to be informative. Use **Previous date**, **Next date** and **Play** to compare available frames; adjust **Speed** to control playback. The date and frame details identify the epoch, approximate wavelength and number of exposures where available. A coadd date represents its grouped observations; inspect the MJD range when precise timing matters.

The **Sky reference map** is a separate optical survey for orientation. It helps locate the field, but its colours, passbands and image resolution differ from the SPHEREx observations.

## Understand overlays and image quality

Under **Overlays on this image**, turn layers on as needed. Some controls are unavailable when their data are absent or in the combined trail view.

- **Show flags:** orange identifies pixels flagged as saturated or contaminated. Flagged regions can affect source identification and measurement.
- **Show expected position:** the green ring marks a catalogue/model-predicted position when a prediction is available. It is not an image measurement. A coordinate-only request describes a fixed field and may have no motion prediction.
- **Show centroid / candidates:** a yellow dashed ring marks an algorithmic centroid candidate; blue dashed circles can show multiple candidates. These are not necessarily the target. Blending, noise, masks and uncertain image alignment can bias the centre. The reported centroid errors do not include all systematic effects.

Violet pixels indicate no usable image data. White fill in masked bright or saturated regions is not a flux measurement. Pixel coverage indicates how much of the field has usable data; even 100% coverage does not establish scientific validation. A pixelated appearance reflects the available image resolution: enlarging the display does not recover finer detail.

Open **Observation details and data quality** for the query centre, source and catalogue assumptions, warnings, frame quality, and **How positions are measured** when available. Read these details before comparing epochs. Depending on the target, catalogue uncertainties, component blending and the adopted coordinate model can limit the prediction.

## Spectrum and light curve

Select **Spectrum** or **Light curve** to request the corresponding graph. Processing can take several minutes. The result may be empty if suitable measurements or comparable spectral responses are unavailable.

- A **spectrum** shows flux measurements versus wavelength. Check the units, measurement quality and wavelength provenance shown with the result. A point at a catalogue position is not, by itself, an independent identification of the target.
- A **light curve** compares measurements over time. Reliable variability claims require suitable photometry and genuinely comparable spectral responses; different wavelengths, bandpasses, masks or calibrations can produce apparent changes.
- **Show preliminary and rejected measurements** reveals diagnostic points when available. These are not equivalent to validated measurements and should not be treated as confirmed detections or connected into a scientific light curve without the required checks.

The project distinguishes validated, provisional and rejected measurement statuses. A validated photometric point is not automatically a validated astrometric position, motion measurement or light curve. Independent reference-star registration and source association are needed for reliable image-based motion measurements.

## Share and revisit a view

Open **Share or embed this view** after choosing a target and settings. Copy the view URL or embed HTML. The URL records the current target or coordinates, field, detector, mode and selected source when applicable; reopening it can require a new archive query. Saved observations refer to a fixed saved selection and are not a fresh archive search.

The [five-target GIF gallery](https://spherex-sky-viewer-zeta.vercel.app/challenge-demo/) presents saved examples and separate assessment results. Its catalogue rings are predictions; the gallery does not establish motion or variability.

## Further reading

- [README](../README.md): project overview, local setup and Python interface.
- [Science remediation and validation](SCIENCE_REMEDIATION.md): validation rules, evidence and limitations.
- [API reference](API.md): request parameters and response fields.
- [Target parameters and sources](../OBJECT_PARAMETERS.md): catalogue origins and assumptions.

Some historical technical reports and saved machine-readable values are in Russian. Data source: NASA SPHEREx via NASA/IPAC IRSA.
