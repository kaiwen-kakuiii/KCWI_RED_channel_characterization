#!/usr/bin/env python
"""Show, from the raw data, why 2023-11-08 grew a spurious 25th slice.

The right edge of one slice steps ~5 px inward partway up the detector.  PypeIt
fits smooth polynomial traces to slit edges, and a step is not smooth, so edge
tracing registered the two levels as two separate edges and paired them into a
4 px "slit" -- the strip between them.  That sliver has no light, so the CONTBARS
frame finds no alignment bars in it and calibration aborts.  Rejecting it leaves
the real slice truncated, and a truncated slice makes the astrometric transform
return NaN, which kills the coadd.

Three views of the same fact:

  left    the dome flat itself, spatial 600-700.  The slice edge is the bright
          to dark transition; it moves left near the top.
  top-r   spatial cuts at two spectral rows, below and above the step.
  bot-r   edge position against spectral row, for the affected slice and its
          neighbours -- the neighbours only drift, this one steps.

Usage:
    python plot_1108_slice_defect.py [-o out.png]
"""
import argparse
import os

import numpy as np
from pypeit.images import rawimage
from pypeit.spectrographs.util import load_spectrograph

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "fits", "by_night", "RH1", "2023-11-08",
                   "KR.20231109.58480.36.fits")

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8880"
CRITICAL, GRID = "#d03b3b", "#e6e5e0"
STEP_ROW = 3480


def edge_at(img, lo, hi, start, plateau_lo):
    """Spatial pixel where the flat falls to half its in-slice level."""
    p = np.median(img[lo:hi], axis=0)
    half = np.median(p[plateau_lo:start]) / 2.0
    idx = np.where(p[start:start + 40] < half)[0]
    return start + idx[0] if len(idx) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out",
                    default=os.path.join(ROOT, "rh1_1108_slice_defect.png"))
    args = ap.parse_args()

    spec = load_spectrograph("keck_kcrm")
    par = spec.default_pypeit_par()["calibrations"]["traceframe"]
    for k in ["use_biasimage", "use_darkimage", "use_pixelflat", "use_illumflat",
              "use_specillum", "mask_cr", "subtract_scattlight"]:
        if k in par["process"].keys():
            par["process"][k] = False
    img = rawimage.RawImage(RAW, spec, 1).process(par["process"]).image
    nspec = img.shape[0]

    fig = plt.figure(figsize=(13.4, 6.4), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.5], hspace=0.42, wspace=0.24)
    axI = fig.add_subplot(gs[:, 0])
    axP = fig.add_subplot(gs[0, 1])
    axE = fig.add_subplot(gs[1, 1])

    # ---- left: the flat itself ------------------------------------------
    x0, x1 = 600, 700
    sub = img[:, x0:x1]
    step = 4                                    # row-average for display only
    disp = sub[: (nspec // step) * step].reshape(-1, step, x1 - x0).mean(axis=1)
    vmin, vmax = np.percentile(disp, [2, 99])
    axI.imshow(disp, origin="lower", aspect="auto", cmap="gray",
               vmin=vmin, vmax=vmax,
               extent=[x0, x1, 0, nspec])
    axI.axhline(STEP_ROW, color=CRITICAL, lw=1.4, ls="--")
    axI.annotate(f"edge steps ~5 px inward\nat spectral row ~{STEP_ROW}",
                 (x1 - 3, STEP_ROW), color=CRITICAL, fontsize=9,
                 ha="right", va="bottom",
                 textcoords="offset points", xytext=(0, 6))
    axI.set_xlabel("spatial pixel", color=INK2, fontsize=10)
    axI.set_ylabel("spectral pixel", color=INK2, fontsize=10)
    axI.set_title("Dome flat, 2023-11-08", color=INK, fontsize=12, loc="left", pad=8)
    axI.tick_params(colors=INK2, labelsize=9)

    # ---- top right: spatial cuts either side of the step -----------------
    for rows, col, lab in [((2000, 2400), MUTED, "row 2000-2400 (below step)"),
                           ((3800, 4200), CRITICAL, "row 3800-4128 (above step)")]:
        lo, hi = rows[0], min(rows[1], nspec)
        p = np.median(img[lo:hi], axis=0)
        axP.plot(np.arange(640, 680), p[640:680], color=col, lw=2, label=lab)
    axP.axvline(653, color=CRITICAL, lw=1, ls=":")
    axP.axvline(658, color=MUTED, lw=1, ls=":")
    axP.annotate("653", (653, axP.get_ylim()[1] * 0.55), color=CRITICAL,
                 fontsize=9, ha="right")
    axP.annotate("658", (658, axP.get_ylim()[1] * 0.72), color=MUTED,
                 fontsize=9, ha="left")
    axP.set_xlabel("spatial pixel", color=INK2, fontsize=10)
    axP.set_ylabel("flat counts", color=INK2, fontsize=10)
    axP.set_title("The edge sits 5 px further in at the red end",
                  color=INK, fontsize=12, loc="left", pad=8)
    axP.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper right")
    axP.tick_params(colors=INK2, labelsize=9)

    # ---- bottom right: edge position vs spectral row ---------------------
    bands = np.arange(0, nspec - 200, 200)
    targets = [("slice 2 right edge", 495, 455, MUTED),
               ("slice 4 right edge", 800, 760, MUTED),
               ("slice 3 right edge (defective)", 640, 600, CRITICAL)]
    for lab, start, pl, col in targets:
        ys = [edge_at(img, b, b + 200, start, pl) for b in bands]
        ys = np.array(ys, dtype=float)
        ys -= np.nanmedian(ys[:6])              # show motion, not absolute position
        axE.plot(bands + 100, ys, "-o", ms=3.5, lw=2, color=col, label=lab,
                 alpha=1.0 if col == CRITICAL else 0.55,
                 zorder=4 if col == CRITICAL else 3)
    axE.axvline(STEP_ROW, color=CRITICAL, lw=1, ls="--", alpha=0.6)
    axE.set_xlabel("spectral pixel", color=INK2, fontsize=10)
    axE.set_ylabel("edge motion (px)", color=INK2, fontsize=10)
    axE.set_title("Neighbours drift; this one steps",
                  color=INK, fontsize=12, loc="left", pad=8)
    axE.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower left")
    axE.tick_params(colors=INK2, labelsize=9)

    for ax in (axP, axE):
        ax.grid(True, color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#d8d7d2")
        ax.set_facecolor(SURFACE)

    fig.savefig(args.out, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
