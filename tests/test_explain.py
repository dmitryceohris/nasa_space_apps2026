"""'Explain this region': verdicts on synthetic epoch grids, the approved wording, and the HTTP path."""
import numpy as np
import pytest

from api import explain_text
from processing.explain import analyse
from processing.frames import FLAG_CODE_NAMES

N_PX, SCALE = 41, 6.15
MJD = np.array([60800.0, 60980.0, 61165.0, 61345.0])  # ~ every six months
GHOST = 1 << FLAG_CODE_NAMES.index("GHOST")
OVERFLOW = 1 << FLAG_CODE_NAMES.index("OVERFLOW")


def star(img, x, y, peak, sigma=1.0):
    yy, xx = np.indices(img.shape)
    img += peak * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma ** 2))


def grids(sources, n=4, noise=1.0, seed=0):
    """sources: callables (epoch index, image, codes) drawing into one epoch."""
    rng = np.random.default_rng(seed)
    img = rng.normal(0, noise, (n, N_PX, N_PX))
    codes = np.zeros((n, N_PX, N_PX), np.uint16)
    for j in range(n):
        for s in sources:
            s(j, img[j], codes[j])
    return {"image": img, "variance": np.full(img.shape, noise ** 2), "codes": codes, "mjd": MJD[:n],
            "wave_min": np.full(n, 2.0), "wave_max": np.full(n, 2.01), "channel_width": np.array(0.05),
            "predicted": np.full((n, 2), np.nan)}


BOX = (10, 10, 30, 30)


def verdicts(facts):
    return [it["verdict"] for it in facts["items"]]


def test_a_star_moving_north_is_a_moving_candidate():
    years = (MJD - MJD[0]) / 365.25
    a = grids([lambda j, im, c: star(im, 20.0, 18.0 + 2.0 * years[j], 60)])
    facts = analyse(a, BOX, SCALE)
    assert verdicts(facts) == ["moving"]
    it = facts["items"][0]
    assert it["direction"] == "north"
    assert it["pm_arcsec_yr"] == pytest.approx(2.0 * SCALE, rel=0.2)
    assert it["n_pair"] >= 2 and it["prediction"] == "none" and it["fast"]


def test_moving_star_on_its_predicted_path_agrees():
    years = (MJD - MJD[0]) / 365.25
    a = grids([lambda j, im, c: star(im, 20.0, 18.0 + 2.0 * years[j], 60)])
    a["predicted"] = np.column_stack([np.full(4, 20.0), 18.0 + 2.0 * years])
    facts = analyse(a, BOX, SCALE, target="Test Star", parallax_mas=10.0)
    it = facts["items"][0]
    assert it["prediction"] == "agrees" and it["offset_px"] < 1.0 and it["target"] == "Test Star"


def test_a_star_that_changes_brightness_in_place():
    a = grids([lambda j, im, c: star(im, 20, 20, 60 if j % 2 else 20)])
    facts = analyse(a, BOX, SCALE)
    assert verdicts(facts) == ["brightness"]
    assert facts["items"][0]["change_pct"] == pytest.approx(100 * 40 / 60, abs=8)


def test_a_constant_star_gives_nothing_significant():
    facts = analyse(grids([lambda j, im, c: star(im, 20, 20, 60)]), BOX, SCALE)
    assert verdicts(facts) == ["nothing"]


def test_a_source_in_one_epoch_only_is_seen_once():
    a = grids([lambda j, im, c: star(im, 15, 25, 60) if j == 2 else None])
    facts = analyse(a, BOX, SCALE)
    assert verdicts(facts) == ["seen_once"]
    assert facts["items"][0]["mjd_seen"] == MJD[2] and facts["items"][0]["n_other"] == 3


def test_a_one_epoch_source_on_ghost_pixels_is_an_artefact():
    def ghost(j, im, c):
        if j == 2:
            star(im, 15, 25, 60)
            c[20:31, 10:21] |= GHOST
    facts = analyse(grids([ghost]), BOX, SCALE)
    assert verdicts(facts) == ["artefact"] and facts["items"][0]["flag_names"] == ["GHOST"]


def test_a_saturated_star_cannot_be_measured():
    def saturated(j, im, c):
        star(im, 20, 20, 200 if j % 2 else 100)
        im[19:22, 19:22] = np.nan
        c[19:22, 19:22] |= OVERFLOW
    facts = analyse(grids([saturated]), BOX, SCALE)
    assert facts["items"][0]["verdict"] == "cant_tell" and facts["items"][0]["reason"] == "saturated"


def test_an_area_without_data_in_two_epochs():
    a = grids([])
    a["image"][1:] = np.nan
    facts = analyse(a, BOX, SCALE)
    assert verdicts(facts) == ["not_enough_data"] and facts["items"][0]["n_epochs_with_data"] == 1


def test_same_region_gives_the_same_answer():
    years = (MJD - MJD[0]) / 365.25
    a = grids([lambda j, im, c: star(im, 20.0, 18.0 + 2.0 * years[j], 60)])
    assert explain_text.render(analyse(a, BOX, SCALE)) == explain_text.render(analyse(a, BOX, SCALE))


# ------------------------------------------------------------------ wording (docs/EXPLAIN_TEXTS.md)
def card(item, n_partial=0, n_with=4):
    return explain_text.render({"n_epochs": 4, "n_epochs_with_data": n_with, "n_partial": n_partial,
                                "items": [item]})["cards"][0]


def test_moving_card_has_the_approved_text():
    c = card({"verdict": "moving", "n_pair": 4, "n_epochs": 5, "mjd_first": 60902.8654, "mjd_last": 61268.6315,
              "shift_px": 1.72, "shift_arcsec": 10.6, "pm_arcsec_yr": 10.33, "pm_err_arcsec_yr": 0.81,
              "direction": "north", "prediction": "agrees", "offset_px": 0.42, "target": "Barnard's Star",
              "small_shift": True, "n_flagged": 0})
    assert c["title"] == "Moving object candidate"
    assert c["paragraphs"] == [
        "Something here is on the move. In 4 of 5 [epochs|epoch] a bright spot appears in a new place and leaves "
        "a dark gap where it used to be.",
        "Between Aug 15, 2025 and Aug 16, 2026 it shifted by 10.6″ (1.7 [pixels|pixel]), about 10.3 ± 0.8″ per "
        "year toward the north.",
        "This agrees with the predicted path of Barnard's Star to within 0.4 px."]
    assert c["caveats"] == [
        "The shift is smaller than two pixels, so the bright and dark spots partly overlap: treat the speed as "
        "approximate.",
        "Image alignment was taken from the telescope's coordinates and was not checked against reference stars."]
    assert c["next"] == "Open the Trail view: the motion shows as a blue-to-yellow streak."


def test_singular_forms_and_flag_names():
    c = card({"verdict": "not_enough_data", "n_epochs_with_data": 1})
    assert c["paragraphs"] == ["This area has usable data in only 1 [epoch]. At least two are needed to see any "
                               "change."]
    c = card({"verdict": "artefact", "flag_names": ["GHOST", "GHOST_FPA", "CROSSTALK"]}, n_partial=1)
    assert "carry the [ghost] and [crosstalk] [flags]" in c["paragraphs"][0]
    assert c["caveats"] == ["Part of this area has no data in 1 epoch."]


def test_brightness_card_rounds_to_five_percent_and_warns_about_wavelength():
    c = card({"verdict": "brightness", "change_pct": 33.4, "mjd_bright": 60800.0, "mjd_faint": 60980.0,
              "n_flagged": 2, "different_wavelengths": True, "wave_min": 1.63, "wave_max": 2.4})
    assert c["paragraphs"][0] == ("The source stays in place, but its brightness changes by about 35%: brightest "
                                  "on May 5, 2025, faintest on Nov 1, 2025.")
    assert c["caveats"][0].startswith("The epochs were observed at different wavelengths (1.63–2.40 µm).")
    assert c["caveats"][1] == "In 2 epochs the source has flagged pixels; those epochs were not used for the numbers above."


def test_header_only_with_several_cards_and_footer_always():
    out = explain_text.render({"n_epochs": 4, "n_epochs_with_data": 4, "n_partial": 0,
                               "items": [{"verdict": "nothing"}]})
    assert out["header"] is None
    assert out["footer"] == "This is an automatic reading of the data, not a confirmed discovery."
    out = explain_text.render({"n_epochs": 4, "n_epochs_with_data": 4, "n_partial": 0,
                               "items": [{"verdict": "cant_tell", "reason": "edge"}] * 2})
    assert out["header"] == "We found 2 sources in this area. Each one is explained separately below."
    assert all(t in explain_text.GLOSSARY for t in ("epoch", "flags", "saturated", "brown dwarfs"))


def test_a_saturated_star_with_bloom_flags_is_not_called_an_artefact():
    bloom = 1 << FLAG_CODE_NAMES.index("BLOOM")

    def star_with_bloom(j, im, c):  # Barnard's Star in real data: NaN core, OVERFLOW, BLOOM all around
        star(im, 20, 20, 200 if j % 2 else 100)
        c[16:25, 16:25] |= bloom
        im[19:22, 19:22] = np.nan
        c[19:22, 19:22] |= OVERFLOW
    facts = analyse(grids([star_with_bloom]), BOX, SCALE)
    assert [(it["verdict"], it.get("reason")) for it in facts["items"]] == [("cant_tell", "saturated")]
    assert facts["n_partial"] == 0  # the saturated core is not "missing data"
