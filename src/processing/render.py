"""numpy -> PNG/GIF.

Соглашение об осях (важно для маркеров на фронтенде):
  в FITS/numpy ось y направлена ВВЕРХ (строка 0 = юг), в PNG — ВНИЗ.
  Поэтому массив всегда переворачивается (np.flipud) и пиксельные координаты для JSON
  пересчитываются через fits_to_png_xy. Оба действия — только здесь.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image


def fits_to_png_xy(x: float, y: float, height: int) -> list[float]:
    return [round(float(x), 3), round(float(height - 1 - y), 3)]


def shared_limits(images, lo: float = 1.0, hi: float = 99.7) -> tuple[float, float]:
    """ОДНИ пределы растяжки на все кадры. Если считать их для каждого кадра отдельно,
    яркость будет «прыгать» и пользователь увидит ложную переменность."""
    vals = np.concatenate([im[np.isfinite(im)].ravel() for im in images if np.isfinite(im).any()] or [np.zeros(1)])
    vmin, vmax = np.percentile(vals, [lo, hi])
    if vmax <= vmin:
        vmax = vmin + 1e-6
    return float(vmin), float(vmax)


def _asinh(img: np.ndarray, vmin: float, vmax: float, soft: float = 0.1) -> np.ndarray:
    x = np.clip((np.nan_to_num(img, nan=vmin) - vmin) / (vmax - vmin), 0, 1)
    return np.arcsinh(x / soft) / np.arcsinh(1 / soft)


def _save(arr_hwc_or_hw: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.flipud(arr_hwc_or_hw)).save(path, optimize=True)


NODATA_RGB = (72, 36, 110)  # dark violet: pixel without data (outside coverage or masked)


def save_gray(img: np.ndarray, path: Path, vmin: float, vmax: float,
              saturated: np.ndarray | None = None) -> None:
    g = (_asinh(img, vmin, vmax) * 255).astype(np.uint8)
    rgb = np.repeat(g[..., None], 3, axis=2)  # the sky stays grey: R = G = B
    nodata = ~np.isfinite(img)
    # no data: a colour no sky pixel can have. The former grey level 25 lay inside the noise of a dark sky
    # (asinh stretch puts the sky at ~0-20), so a gap in coverage looked like empty sky.
    rgb[nodata] = NODATA_RGB
    if saturated is not None:
        rgb[nodata & saturated] = 255  # насыщенное ядро яркой звезды — белым, а не «дырой»
    _save(rgb, path)


def save_diff(sig: np.ndarray, path: Path, lim: float = 6.0) -> None:
    """Дивергентная палитра: синий = стало темнее, красный = стало ярче, белый = без изменений."""
    from matplotlib import colormaps

    x = np.clip(np.nan_to_num(sig) / lim, -1, 1) * 0.5 + 0.5
    rgb = (colormaps["RdBu_r"](x)[..., :3] * 255).astype(np.uint8)
    rgb[~np.isfinite(sig)] = (40, 40, 40)
    _save(rgb, path)


def save_rgb(rgb01: np.ndarray, path: Path, background: np.ndarray | None = None,
             vmin: float = 0, vmax: float = 1) -> None:
    """«След»: цветные остатки поверх серого шаблона (шаблон приглушён до 35%)."""
    base = np.zeros_like(rgb01)
    if background is not None:
        base = 0.35 * _asinh(background, vmin, vmax)[..., None] * np.ones(3)
    out = np.clip(base + rgb01, 0, 1)
    _save((out * 255).astype(np.uint8), path)


def save_mask(mask: np.ndarray, path: Path, color=(255, 80, 0), alpha: int = 140) -> None:
    """Прозрачный PNG-слой: отмеченные пиксели — полупрозрачным оранжевым."""
    rgba = np.zeros((*mask.shape, 4), np.uint8)
    rgba[mask] = (*color, alpha)
    _save(rgba, path)


def save_gif(png_paths: list[Path], out: Path, ms_per_frame: int = 600, scale: int = 8) -> None:
    frames = [Image.open(p).convert("RGB") for p in png_paths]
    frames = [f.resize((f.width * scale, f.height * scale), Image.NEAREST) for f in frames]
    out.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=ms_per_frame, loop=0)
