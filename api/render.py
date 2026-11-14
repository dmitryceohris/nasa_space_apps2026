"""GIF export and the PNG pixel convention. Frames and masks are drawn by src/processing/render.py.

Convention (same as there): PNG has y down, FITS/numpy has y up; every pixel coordinate sent to the front goes
through png_xy() / processing.render.fits_to_png_xy(), otherwise markers drift off the objects."""
from PIL import Image


def png_xy(x, y, height):
    """FITS-orientation pixel (0-based, y up) -> PNG pixel (0-based, y down)."""
    return [float(x), float(height - 1 - y)]


def save_gif(paths, out, ms=600):
    """`out` may be a file or a BytesIO."""
    ims = [Image.open(p).convert("P") for p in paths]
    ims[0].save(out, format="GIF", save_all=True, append_images=ims[1:], duration=ms, loop=0)
