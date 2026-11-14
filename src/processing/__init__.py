"""Обработка снимков SPHEREx: эпохи, общая сетка, coadd/diff/trail, PNG, спектр."""
from .epochs import group_epochs, select_per_epoch
from .frames import Frame, coadd, difference, make_grid, reproject_exposure, trail
from .photometry import forced_photometry, light_curve, spectrum
from .pipeline import (build_spectrum, build_view, frames_from_exposures,
                       lightcurve_from_spectrum, render_view)
from .targets import PRESETS, Target, find_preset, load_targets, resolve_simbad

__all__ = [
    "group_epochs", "select_per_epoch", "Frame", "coadd", "difference", "make_grid",
    "reproject_exposure", "trail", "forced_photometry", "light_curve", "spectrum",
    "build_view", "build_spectrum", "frames_from_exposures", "render_view",
    "lightcurve_from_spectrum", "PRESETS", "Target", "find_preset", "load_targets", "resolve_simbad",
]
