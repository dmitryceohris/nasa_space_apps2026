"""Approved wording of the "Explain this region" cards (docs/EXPLAIN_TEXTS.md). One fixed text per block.

Markup: [term] or [shown text|glossary key] is an underlined term with a tooltip from GLOSSARY.
"""
from astropy.time import Time

GLOSSARY = {
    "epoch": "A group of SPHEREx images of the same place taken within a few days. SPHEREx returns to each part "
             "of the sky about every six months.",
    "pixel": "The smallest element of the image. One SPHEREx pixel covers 6.15″ of sky.",
    "arcsecond": "A unit of angle: 1/3600 of a degree. The full Moon is about 1800″ across.",
    "proper motion": "How fast a star drifts across the sky, in arcseconds per year. Nearby stars drift fastest.",
    "parallax": "A small yearly back-and-forth shift of nearby objects, caused by the Earth's orbit around the Sun.",
    "residual": "What is left after subtracting the other epochs. Bright means more light now, dark means less.",
    "saturated": "So bright that the detector reached its limit, so the true brightness and centre are unknown.",
    "flags": "Labels that SPHEREx attaches to pixels affected by instrument problems.",
    "ghost": "A faint false image of a bright star, caused by reflections inside the telescope.",
    "crosstalk": "A false signal copied from a bright source into another part of the detector by the electronics.",
    "cosmic ray": "A charged particle that hits the detector and lights up one or a few pixels for a single image.",
    "noise": "Random variation in every image. Changes smaller than the noise cannot be told apart from chance.",
    "S/N": "Signal-to-noise ratio: how many times a signal is stronger than the noise. Below 5 we do not trust it.",
    "spectrum": "How bright an object is at each wavelength of light.",
    "wavelength": "The \"colour\" of light. SPHEREx sees infrared light from 0.75 to 5.0 µm, invisible to the eye.",
    "variability": "A real change in the object's own brightness over time.",
    "variable stars": "Stars whose brightness changes, for example pulsating Mira stars or young stars still "
                      "gathering gas.",
    "brown dwarfs": "Objects heavier than planets but too light to shine like stars. They are cold and glow mainly "
                    "in infrared light.",
    "asteroids": "Small rocky bodies orbiting the Sun. Seen from Earth they move very fast across the sky.",
    "detector": "SPHEREx has six detectors; each covers its own range of infrared wavelengths.",
}

FLAG_SHOWN = {"GHOST": "ghost", "GHOST_FPA": "ghost", "GHOST_EXT": "ghost", "CROSSTALK": "crosstalk",
              "STREAK": "streak", "BLOOM": "bloom", "SNOWBALL": "snowball", "HALO": "halo",
              "OVERFLOW": "saturation", "NONLINEAR": "nonlinear response"}

HEADER = "We found {n_sources} sources in this area. Each one is explained separately below."
FOOTER = "This is an automatic reading of the data, not a confirmed discovery."

TITLES = {"not_enough_data": "Not enough data", "artefact": "Likely artefact", "cant_tell": "Can't tell",
          "moving": "Moving object candidate", "seen_once": "Seen only once", "brightness": "Brightness change",
          "nothing": "Nothing significant"}

NEXT = {
    "not_enough_data": "Try a larger field of view or another [detector].",
    "artefact": "Turn on the flags layer to see the affected pixels.",
    "cant_tell": "Look at the Co-add view: the source is real, only the measurement is limited.",
    "moving": "Open the Trail view: the motion shows as a blue-to-yellow streak.",
    "seen_once": "Check this date in the Single view: a real object has the shape of a star, a cosmic ray is one "
                 "sharp pixel.",
    "brightness": "Open the Light curve tab to follow the brightness at one wavelength.",
    "nothing": "Select a smaller area around a source, or try the Difference view.",
}

CANT_TELL = {
    "saturated": "The centre of this source is [saturated]: it is so bright that the detector reached its limit, so "
                 "its exact position and brightness cannot be measured.",
    "ambiguous": "Two sources are too close to tell which one changes.",
    "edge": "The source touches the edge of the data, so its centre would be biased.",
}


def plural(n, one, many=None):
    """'1 epoch' / '2 epochs'. one/many may carry glossary markup."""
    return f"{n} {one if n == 1 else (many or one + 's')}"


def date(mjd: float) -> str:
    d = Time(mjd, format="mjd").to_datetime()
    return f"{d:%b} {d.day}, {d.year}"


def num(v: float) -> str:
    return f"{v:.1f}"


def flag_names(names) -> str:
    shown = []
    for n in names:
        s = FLAG_SHOWN.get(n, n.lower())
        s = f"[{s}]" if s in GLOSSARY else s
        if s not in shown:
            shown.append(s)
    if len(shown) <= 1:
        return "".join(shown)
    return ", ".join(shown[:-1]) + " and " + shown[-1]


def _paragraphs(it: dict) -> list[str]:
    v = it["verdict"]
    if v == "not_enough_data":
        return [f"This area has usable data in only {plural(it['n_epochs_with_data'], '[epoch]', '[epochs|epoch]')}. "
                "At least two are needed to see any change."]
    if v == "artefact":
        return [f"This is probably not a real object. The pixels here carry the {flag_names(it['flag_names'])} "
                "[flags], which SPHEREx uses to mark instrument effects rather than light from the sky."]
    if v == "cant_tell":
        return [CANT_TELL[it["reason"]]]
    if v == "moving":
        shift = num(it["shift_px"])
        pixels = f"{shift} " + ("[pixel]" if shift == "1.0" else "[pixels|pixel]")
        out = [f"Something here is on the move. In {it['n_pair']} of {it['n_epochs']} [epochs|epoch] a bright spot "
               "appears in a new place and leaves a dark gap where it used to be.",
               f"Between {date(it['mjd_first'])} and {date(it['mjd_last'])} it shifted by {num(it['shift_arcsec'])}″ "
               f"({pixels}), about {num(it['pm_arcsec_yr'])} ± {num(it['pm_err_arcsec_yr'])}″ per year toward "
               f"the {it['direction']}."]
        p = it.get("prediction", "none")
        if p == "agrees":
            out.append(f"This agrees with the predicted path of {it['target']} to within {num(it['offset_px'])} px.")
        elif p == "disagrees":
            out.append(f"The shift does not match the predicted path of {it['target']} "
                       f"(off by {num(it['offset_px'])} px). This may be a different object.")
        elif it.get("fast"):
            out.append("Objects that drift this fast are usually nearby stars or [brown dwarfs]: the closer an "
                       "object is, the faster it seems to move. This is a measured [proper motion] candidate, "
                       "not a confirmed one.")
        else:
            out.append("The shift is small. It may be real motion, or the images may be lined up slightly "
                       "differently.")
        return out
    if v == "seen_once":
        return [f"A source appears here on {date(it['mjd_seen'])} and is absent in the other "
                f"{plural(it['n_other'], '[epoch]', '[epochs|epoch]')}.",
                "Fast movers such as [asteroids] and comets cross this field in days, so they show up in just one "
                "epoch. A [cosmic ray] or an instrument effect can look the same."]
    if v == "brightness":
        pct = max(5, int(round(it["change_pct"] / 5)) * 5)
        return [f"The source stays in place, but its brightness changes by about {pct}%: brightest on "
                f"{date(it['mjd_bright'])}, faintest on {date(it['mjd_faint'])}.",
                "This is how [variable stars] look in the data."]
    if v == "nothing":
        return [f"Nothing in this area changes beyond the [noise] level across "
                f"{plural(it['n_epochs'], '[epoch]', '[epochs|epoch]')}."]
    raise ValueError(v)


def _caveats(it: dict, n_partial: int) -> list[str]:
    v, out = it["verdict"], []
    if v == "not_enough_data":
        return out
    if v in ("moving", "brightness") and it.get("different_wavelengths"):
        out.append(f"The epochs were observed at different wavelengths ({it['wave_min']:.2f}–{it['wave_max']:.2f} µm). "
                   "Part of the change may come from the object's [spectrum], not from real [variability].")
    if v == "moving" and it.get("small_shift"):
        out.append("The shift is smaller than two pixels, so the bright and dark spots partly overlap: treat the "
                   "speed as approximate.")
    if v == "moving" and it.get("parallax_noticeable"):
        out.append("Part of this shift is [parallax]: a yearly wobble caused by the Earth moving around the Sun.")
    if v in ("moving", "brightness") and it.get("n_flagged"):
        out.append(f"In {plural(it['n_flagged'], 'epoch')} the source has flagged pixels; those epochs were not "
                   "used for the numbers above.")
    if n_partial:
        out.append(f"Part of this area has no data in {plural(n_partial, 'epoch')}.")
    if v == "moving":
        out.append("Image alignment was taken from the telescope's coordinates and was not checked against "
                   "reference stars.")
    return out


def render(facts: dict) -> dict:
    """facts from processing.explain.analyse -> {header, cards: [{verdict, title, paragraphs, caveats, next, xy}],
    footer, glossary}."""
    cards = []
    for it in facts["items"]:
        it = dict(it)
        if it["verdict"] == "nothing":
            it["n_epochs"] = facts["n_epochs_with_data"]
        cards.append({"verdict": it["verdict"], "title": TITLES[it["verdict"]], "paragraphs": _paragraphs(it),
                      "caveats": _caveats(it, facts["n_partial"]), "next": NEXT[it["verdict"]],
                      "xy": it.get("xy")})
    return {"header": HEADER.format(n_sources=len(cards)) if len(cards) > 1 else None, "cards": cards,
            "footer": FOOTER, "glossary": GLOSSARY, "facts": facts}
