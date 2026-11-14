"""Quick look: a blinking animation of a Stack (in Jupyter it plays inline)."""
from __future__ import annotations

import numpy as np

from .stack import Stack


def show(st: Stack, mode: str = "coadd", interval: int = 600, track=None, predicted: bool = True,
         size: float = 5.0):
    """Animation of the epochs. mode: 'coadd' (shared asinh stretch, so brightness does not jump) or 'diff'
    (red = brighter than the other epochs, blue = fainter). North up, east left.
    track: a Track from spx.track() to mark the measured centre (yellow ring); predicted: green ring at the
    catalogue position. In Jupyter returns an inline HTML player, elsewhere a matplotlib FuncAnimation
    (call .save('movie.gif') or plt.show())."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    from processing.render import shared_limits

    from .measure import predict
    from .stack import diff as _diff

    if not len(st):
        raise ValueError("пустой стек")
    if mode == "diff":
        cube, cmap, vmin, vmax = _diff(st), "RdBu_r", -6.0, 6.0
    elif mode == "coadd":
        lo, hi = shared_limits(list(st.cube))
        soft = 0.1
        cube = np.arcsinh(np.clip((st.cube - lo) / (hi - lo), 0, 1) / soft) / np.arcsinh(1 / soft)
        cmap, vmin, vmax = "gray", 0.0, 1.0
    else:
        raise ValueError("mode: 'coadd' или 'diff'")

    fig, ax = plt.subplots(figsize=(size, size))
    im = ax.imshow(np.ma.masked_invalid(cube[0]), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    im.cmap.set_bad("#282828")
    ax.set_xlabel("x, px (east ←)")
    ax.set_ylabel("y, px (north ↑)")
    has_model = bool(st.target.pmra or st.target.pmdec or st.target.parallax)
    pred_pts = [st.pixel(predict(st.target, m)) for m in st.mjd] if (predicted and has_model) else None
    pred_mark, = ax.plot([], [], "o", ms=18, mfc="none", mec="#2c2", mew=1.5, label="prediction")
    meas_mark, = ax.plot([], [], "o", ms=12, mfc="none", mec="#fd0", mew=1.5, label="measured")
    if pred_pts or track is not None:
        ax.legend(loc="upper right", fontsize=8, framealpha=0.6)
    title = ax.set_title("")

    def draw(i):
        im.set_data(np.ma.masked_invalid(cube[i]))
        if pred_pts:
            pred_mark.set_data([pred_pts[i][0]], [pred_pts[i][1]])
        if track is not None and track[i]["xy"] is not None:
            meas_mark.set_data([track[i]["xy"][0]], [track[i]["xy"][1]])
        else:
            meas_mark.set_data([], [])
        wl = f", {st.wave[i]:.2f} µm" if np.isfinite(st.wave[i]) else ""
        title.set_text(f"{st.target.name}: epoch {i + 1}/{len(st)}  {st.dates[i]}{wl}  [{mode}]")
        return im, pred_mark, meas_mark, title

    draw(0)
    anim = FuncAnimation(fig, draw, frames=len(st), interval=interval, blit=False)
    try:
        from IPython import get_ipython
        from IPython.display import HTML

        if get_ipython() is not None and "IPKernelApp" in get_ipython().config:
            html = HTML(anim.to_jshtml())
            plt.close(fig)
            return html
    except ImportError:
        pass
    return anim
