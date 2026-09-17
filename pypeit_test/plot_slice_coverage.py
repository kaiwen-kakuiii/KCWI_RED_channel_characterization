#!/usr/bin/env python
"""Why the ends of a KCRM datacube get trimmed: the 24 slices do not agree on λ.

Each slice of the image slicer meets the grating at a slightly different angle,
so each records a slightly different wavelength window.  The datacube spans the
UNION of those windows, which means that near its ends only SOME slices carry
data -- and the coadd fills the missing spaxels by smearing flux in from the
slices that do, because weighted-mean resampling is only flux-conserving where
coverage is complete.  The counts there look completely normal.

The figure is the argument: one bar per slice, and the interval every slice
covers marked on top of them.  Everything outside that interval is a wavelength
where the cube is partly fabricated.

Windows are measured the way measure_trim.py measures them -- per slit, the
min/max good-pixel wavelength of WAVEIMG inside the slit edges, BPM applied.

Usage:
    python plot_slice_coverage.py [--spec2d <file>] [--outdir <dir>]
"""
import argparse
import glob
import os

import numpy as np
from astropy.io import fits

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402

import plot_blaze_summary as P

ROOT = P.ROOT
INK, INK2 = P.INK, P.INK2
GRID, SPINE, SURFACE = P.GRID, P.SPINE, P.SURFACE
BAR = P.HUES[0]
WARN = P.HUES[1]

DEFAULT = os.path.join(
    ROOT, "RL pypeit run", "2024-12-28", "pypeit_run", "keck_kcrm_B",
    "Science", "spec2d_*.fits")


def windows(path):
    """(lo, hi) per slit, in Angstrom, exactly as measure_trim.py reads them."""
    with fits.open(path, memmap=True) as h:
        pick = lambda tag: next(x.data for x in h if x.name.endswith(tag))
        wave, bpm, slits = pick("WAVEIMG"), pick("BPMMASK"), pick("SLITS")
        xx = np.arange(wave.shape[1])[None, :]
        out = []
        for i in range(len(slits)):
            l = int(np.median(slits["left_init"][i]))
            r = int(np.median(slits["right_init"][i]))
            m = (xx >= l) & (xx <= r) & (wave > 0) & (bpm == 0)
            wv = wave[m]
            out.append((float(wv.min()), float(wv.max())))
    return out


def figure(win, path, titles=False):
    lo = max(a for a, _ in win)          # full coverage starts at the LAST start
    hi = min(b for _, b in win)          # and ends at the FIRST end
    umin = min(a for a, _ in win)
    umax = max(b for _, b in win)
    order = np.argsort([a for a, _ in win])

    fig, ax = plt.subplots(figsize=(13.33, 5.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    n = len(win)
    # The partial zones: wavelengths the cube covers but not every slice does.
    ax.axvspan(umin, lo, color=WARN, alpha=0.13, lw=0, zorder=1)
    ax.axvspan(hi, umax, color=WARN, alpha=0.13, lw=0, zorder=1)
    ax.axvspan(lo, hi, color=BAR, alpha=0.07, lw=0, zorder=1)

    for y, i in enumerate(order):
        a, b = win[i]
        ax.plot([a, b], [y, y], color=BAR, lw=4.0, solid_capstyle="butt",
                zorder=3)

    for x in (lo, hi):
        ax.axvline(x, color=INK, lw=1.6, ls=(0, (5, 3)), zorder=4)

    ax.set_ylim(-1.6, n + 4.4)
    ax.set_xlim(umin - 0.045 * (umax - umin), umax + 0.045 * (umax - umin))
    ax.set_yticks([])
    ax.set_ylabel(f"the {n} slices\nsorted by where each starts", fontsize=12,
                  color=INK, linespacing=1.5)
    ax.tick_params(axis="x", colors=INK, labelsize=12, length=4, color=SPINE)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("wavelength (Å)", fontsize=13, color=INK)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(SPINE)

    ax.annotate("every slice covers this", xy=(0.5 * (lo + hi), n + 1.4),
                ha="center", va="center", fontsize=13, color=INK,
                fontweight="bold")
    ax.annotate("", xy=(lo, n + 0.35), xytext=(hi, n + 0.35),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.6))
    ax.text(0.5 * (lo + hi), n + 2.9, f"{hi - lo:.0f} Å kept",
            ha="center", va="center", fontsize=12, color=INK)

    for x0, x1, lab in ((umin, lo, "only SOME slices\nreach here"),
                        (hi, umax, "only SOME slices\nreach here")):
        ax.text(0.5 * (x0 + x1), n + 1.6, lab, ha="center", va="center",
                fontsize=11.5, color=WARN, fontweight="bold", linespacing=1.5)
        ax.text(0.5 * (x0 + x1), -0.9, f"{x1 - x0:.0f} Å", ha="center",
                va="center", fontsize=12, color=WARN, fontweight="bold")

    if titles:
        ax.set_title("The 24 slices do not agree on wavelength", color=INK,
                     fontsize=16, loc="left", fontweight="bold", pad=42)

    fig.subplots_adjust(left=0.115, right=0.985, top=0.93, bottom=0.115)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return lo, hi, umin, umax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec2d", default=DEFAULT)
    ap.add_argument("--outdir", default=ROOT)
    a = ap.parse_args()
    paths = sorted(glob.glob(a.spec2d))
    if not paths:
        raise SystemExit(f"no spec2d matched {a.spec2d}")
    win = windows(paths[0])
    out = os.path.join(a.outdir, "kcrm_slice_coverage.png")
    lo, hi, umin, umax = figure(win, out)
    print(f"{os.path.basename(paths[0])}\n{len(win)} slices, "
          f"union {umin:.0f}-{umax:.0f} A ({umax - umin:.0f} A), "
          f"full coverage {lo:.0f}-{hi:.0f} A ({hi - lo:.0f} A)\n"
          f"partial: blue {lo - umin:.0f} A, red {umax - hi:.0f} A -> {out}")


if __name__ == "__main__":
    main()
