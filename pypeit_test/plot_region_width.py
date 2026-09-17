#!/usr/bin/env python
"""Illustrate how find_object_regions.py sets the star band, on a real frame.

One slice, one profile.  The band is a FLUX-ENCLOSURE interval, not a fitted
profile: clip the sky-subtracted along-slice profile at zero, take the 2% and
98% points of its cumulative, pad 8 percentage points each side.  No Gaussian is
assumed, so an asymmetric or slightly saturated star still gives a sane band.

What is emitted is the COMPLEMENT -- `user_regions = :a,z:` declares the sky, by
excluding the star.

Usage:
    python plot_region_width.py [--raw <file>] [--calib <dir>] [--outdir <dir>]
"""
import argparse
import glob
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402

from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import rawimage
from pypeit.slittrace import SlitTraceSet

import plot_blaze_summary as P

ROOT = P.ROOT
INK, INK2 = P.INK, P.INK2
GRID, SPINE, SURFACE = P.GRID, P.SPINE, P.SURFACE
STAR, SKY = P.HUES[1], P.HUES[0]        # orange = star band, blue = sky

RUN = os.path.join(ROOT, "RH1 pypeit run", "2023-11-17", "pypeit_run",
                   "keck_kcrm_B")
RAW = os.path.join(ROOT, "fits", "by_night", "RH1", "2023-11-17",
                   "KR.20231118.55188.88.fits")
PAD = 8.0


def profile(raw, calib, spec_name="keck_kcrm"):
    """The brightest slice's sky-subtracted along-slice profile, as % of slice."""
    spec = load_spectrograph(spec_name)
    par = spec.default_pypeit_par()["scienceframe"]
    for k in ("use_biasimage", "use_darkimage", "use_pixelflat", "use_illumflat",
              "use_specillum", "mask_cr", "subtract_scattlight"):
        if k in par["process"].keys():
            par["process"][k] = False
    sci = rawimage.RawImage(raw, spec, 1).process(par["process"]).image
    slits = SlitTraceSet.from_file(
        glob.glob(os.path.join(calib, "Slits_*.fits.gz"))[0])

    best, bestflux = None, -np.inf
    for i in range(slits.nslits):
        lo = int(np.median(slits.left_init[:, i]))
        hi = int(np.median(slits.right_init[:, i]))
        p = np.median(sci[:, lo:hi], axis=0)
        resid = p - np.percentile(p, 25)
        f = np.clip(resid, 0, None).sum()
        if f > bestflux:
            best, bestflux = (i, resid, np.percentile(p, 25)), f
    return best


def figure(idx, resid, base, path):
    n = resid.size
    x = 100.0 * np.arange(n) / n
    c = np.clip(resid, 0, None)
    cum = np.cumsum(c) / c.sum()
    lo = 100.0 * np.searchsorted(cum, 0.02) / n
    hi = 100.0 * np.searchsorted(cum, 0.98) / n
    a, z = max(0.0, lo - PAD), min(100.0, hi + PAD)

    fig, ax = plt.subplots(figsize=(11.5, 4.8), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.axvspan(0, a, color=SKY, alpha=0.14, lw=0, zorder=1)
    ax.axvspan(z, 100, color=SKY, alpha=0.14, lw=0, zorder=1)
    ax.axvspan(lo, hi, color=STAR, alpha=0.20, lw=0, zorder=2)
    for v in (a, z):
        ax.axvline(v, color=STAR, lw=2.0, zorder=4)
    for v in (lo, hi):
        ax.axvline(v, color=STAR, lw=1.4, ls=(0, (4, 3)), zorder=4)

    ax.plot(x, resid, color=INK, lw=2.0, zorder=5)
    ax.axhline(0, color=INK2, lw=1.0, zorder=3)

    top = resid.max()
    ax.annotate("", xy=(lo, top * 1.06), xytext=(hi, top * 1.06),
                arrowprops=dict(arrowstyle="<->", color=STAR, lw=1.6))
    ax.text(0.5 * (lo + hi), top * 1.14, "2%–98% of the flux",
            ha="center", va="bottom", fontsize=12, color=STAR,
            fontweight="bold")
    for v, s in ((a, -1), (z, 1)):
        ax.annotate("", xy=(v, top * 0.62), xytext=(v - s * PAD, top * 0.62),
                    arrowprops=dict(arrowstyle="<->", color=STAR, lw=1.3))
    ax.text(a - PAD / 2, top * 0.70, f"+{PAD:.0f}%", ha="center", va="bottom",
            fontsize=11, color=STAR)
    ax.text(z + PAD / 2, top * 0.70, f"+{PAD:.0f}%", ha="center", va="bottom",
            fontsize=11, color=STAR)

    ax.text(a / 2, top * 0.90, "sky", ha="center", va="center", fontsize=14,
            color=SKY, fontweight="bold")
    ax.text((z + 100) / 2, top * 0.90, "sky", ha="center", va="center",
            fontsize=14, color=SKY, fontweight="bold")

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.06 * top, top * 1.30)
    ax.set_xlabel("position along the slice (%)", fontsize=14, color=INK)
    ax.set_ylabel("counts above sky", fontsize=14, color=INK)
    ax.tick_params(colors=INK, labelsize=12, length=4, color=SPINE)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(SPINE)

    ax.set_title(f"user_regions = :{a:.0f},{z:.0f}:   —   the star band is "
                 f"excluded, everything else is sky",
                 fontsize=14, color=INK, loc="left", pad=12)

    fig.subplots_adjust(left=0.085, right=0.985, top=0.885, bottom=0.145)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    return lo, hi, a, z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=RAW)
    ap.add_argument("--calib", default=os.path.join(RUN, "Calibrations"))
    ap.add_argument("--outdir", default=ROOT)
    a = ap.parse_args()
    idx, resid, base = profile(a.raw, a.calib)
    out = os.path.join(a.outdir, "kcrm_region_width.png")
    lo, hi, aa, zz = figure(idx, resid, base, out)
    print(f"slice {idx}: flux 2-98% = {lo:.0f}-{hi:.0f}% ({hi - lo:.0f}% wide), "
          f"padded to {aa:.0f}-{zz:.0f}%  ->  user_regions = :{aa:.0f},{zz:.0f}:")


if __name__ == "__main__":
    main()
