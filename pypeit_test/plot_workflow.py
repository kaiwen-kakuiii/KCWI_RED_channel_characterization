#!/usr/bin/env python
"""One-slide first-pass workflow figure for the KCRM throughput talk.

Deliberately STOCK PypeIt end to end.  This slide is the first pass -- what the
pipeline does when it works -- and the parameter tuning and the problems that
had to be solved (wavelength templates, the coverage trim) get their own slides
later.  Nothing here should need a caveat.

Three tiers, walked in order:

  the SPINE carries the data products -- the things that exist on disk, which is
  what an audience can hold on to: raw frames, spec2d, datacube, spec1d,
  throughput.
  BELOW each arrow, the operation that makes that transition.
  BELOW that, one panel opening up run_pypeit, because "reduce the data" is the
  step an audience has no picture of: it is really nine calibration products
  built from six kinds of frame, and which frame feeds which step is the
  content.

Frame types are the ones PypeIt actually assigns for keck_kcrm, read off a
.pypeit file, not paraphrased: bias / align / arc / tilt / illumflat,trace /
pixelflat,scattlight / science.

Sizes are set for projection: 16:9, nothing below 9 pt.

Usage:
    python plot_workflow.py [--outdir <dir>] [--pdf]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch   # noqa: E402

import plot_blaze_summary as P

ROOT = P.ROOT
INK, INK2 = P.INK, P.INK2
GRID, SPINE, SURFACE = P.GRID, P.SPINE, P.SURFACE

STEP = P.HUES[0]         # the operations
PROD = "#2b2a28"         # data products, on the spine

# Which frame makes which calibration product, in the order PypeIt builds them.
# The order is the point: you cannot flat-field or wavelength-calibrate until
# you know where the 24 slices fall, so the dome flat's trace comes first.
CALIB = [
    ("bias", "bias", "bias and overscan level"),
    ("dome flat", "trace", "slit edges — where the 24 slices fall"),
    ("contbars", "align", "spatial alignment along each slice"),
    ("arc lamp", "arc", "wavelength solution, λ per pixel"),
    ("arc lamp", "tilt", "tilt of λ across each slice"),
    ("flat lamp", "pixelflat", "pixel-to-pixel response"),
    ("dome flat", "illumflat", "illumination across the field"),
    ("flat lamp", "scattlight", "scattered-light model"),
]
SCIENCE = ("science", "the star", "sky model, then subtract → spec2d")


VS = 1.0        # vertical scale: axes fractions per unit of a 7.5in figure


def box(ax, x, y, w, h, title, sub=None, fc="#ffffff", ec=SPINE, tc=INK,
        fs=13.0, sub_c=INK, sub_fs=10.0, lw=1.5):
    """A rounded box centred on (x, y), laid out from its TOP edge.

    Anchoring the title to the top and hanging the sub-block under it is what
    keeps a multi-line caption from growing back over its own title, which is
    what a centred layout does.
    """
    h = h * VS
    top = y + h / 2
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=0.011",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3))
    if sub:
        ax.text(x, top - 0.030 * VS, title, ha="center", va="center",
                fontsize=fs, color=tc, fontweight="bold", zorder=4)
        ax.text(x, top - 0.056 * VS, sub, ha="center", va="top", fontsize=sub_fs,
                color=sub_c, zorder=4, linespacing=1.6)
    else:
        ax.text(x, y, title, ha="center", va="center", fontsize=fs, color=tc,
                fontweight="bold", zorder=4)


def arrow(ax, p0, p1, color=INK2, lw=1.8):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=15, linewidth=lw,
        color=color, zorder=2, shrinkA=2, shrinkB=2))


def figure(path, detail=True, titles=True, figh=7.5):
    """`detail` opens up run_pypeit; `titles` prints the slide's own heading.

    The talk copy wants neither -- it is dropped into a slide that carries its
    own title -- so it is also drawn short, as a strip, rather than as a full
    page with two thirds of it empty.
    """
    global VS
    VS = 7.5 / figh
    fig, ax = plt.subplots(figsize=(13.33, figh), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    SPINE_Y, OP_Y = (0.858, 0.668) if detail else (0.750, 0.260)

    box(ax, 0.056, SPINE_Y, 0.112, 0.096, "KOA", "public archive",
        fc="#ffffff", ec=INK2, tc=INK, fs=13.5, sub_fs=10.0)

    prods = [
        (0.232, "raw frames", "one night, one setup"),
        (0.402, "spec2d", "24 slices per exposure"),
        (0.572, "datacube", "spaxels × wavelength"),
        (0.742, "spec1d", "counts vs wavelength"),
        (0.917, "throughput", "photons in → recorded"),
    ]
    for x, t, sb in prods:
        box(ax, x, SPINE_Y, 0.145, 0.096, t, sb, fc=PROD, ec=PROD,
            tc="#ffffff", fs=14.0, sub_c="#ffffff", sub_fs=10.0)
    arrow(ax, (0.114, SPINE_Y), (0.157, SPINE_Y), lw=2.0)
    for (x0, _, _), (x1, _, _) in zip(prods, prods[1:]):
        arrow(ax, (x0 + 0.0735, SPINE_Y), (x1 - 0.0735, SPINE_Y), lw=2.0)

    ops = [
        (0.136, 0.162, "download script",
         "> 10 s, star present,\nreal setup"),
        (0.317, 0.162, "run_pypeit",
         "calibrate, then\nsky-subtract each frame"),
        (0.487, 0.162, "coadd_datacube",
         "resample all 24 slices\nonto one grid"),
        (0.657, 0.162, "extract_datacube",
         "sum the light around\nthe star, at every λ"),
        (0.8295, 0.162, "sensfunc",
         "counts ÷ standard\n→ zeropoint → fit"),
    ]
    for x, w, t, sb in ops:
        box(ax, x, OP_Y, w, 0.118, t, sb, ec=STEP, tc=STEP, fs=12.5, lw=1.8)
        arrow(ax, (x, OP_Y + 0.059 * VS), (x, SPINE_Y - 0.050 * VS),
              color=STEP, lw=1.6)
    # ---- inside run_pypeit --------------------------------------------------
    if not detail:
        _finish(fig, ax, titles, path)
        return
    PX, PW = 0.500, 0.600
    PTOP, PBOT = 0.575, 0.035
    ax.add_patch(FancyBboxPatch(
        (PX - PW / 2, PBOT), PW, PTOP - PBOT,
        boxstyle="round,pad=0,rounding_size=0.011",
        facecolor="#ffffff", edgecolor=STEP, linewidth=1.8, zorder=1))
    arrow(ax, (0.317, PTOP), (0.317, OP_Y - 0.059), color=STEP, lw=1.6)

    ax.text(PX, PTOP - 0.034, "inside run_pypeit — which frames make which "
            "calibration", ha="center", va="center", fontsize=12.5,
            color=STEP, fontweight="bold", zorder=4)
    ax.text(PX, PTOP - 0.070, "order matters: nothing can be flat-fielded or "
            "wavelength-calibrated\nuntil the slit edges say where the slices "
            "fall", ha="center", va="center", fontsize=10.5, color=INK,
            style="italic", linespacing=1.5, zorder=4)

    LX, TX, RX = 0.245, 0.370, 0.487
    ax.text(LX, PTOP - 0.126, "FRAME", fontsize=10.0, color=INK,
            fontweight="bold", va="center", zorder=4)
    ax.text(RX, PTOP - 0.126, "WHAT IT CALIBRATES", fontsize=10.0, color=INK,
            fontweight="bold", va="center", zorder=4)

    y = PTOP - 0.160
    for name, tag, makes in CALIB:
        ax.text(LX, y, name, fontsize=11.5, color=INK, fontweight="bold",
                va="center", zorder=4)
        ax.text(TX, y, tag, fontsize=10.5, color=INK, va="center",
                fontfamily="monospace", zorder=4)
        ax.text(RX - 0.025, y, "→", fontsize=11.5, color=STEP, va="center",
                ha="center", zorder=4)
        ax.text(RX, y, makes, fontsize=11.5, color=INK, va="center", zorder=4)
        y -= 0.0415

    ax.plot([LX, PX + PW / 2 - 0.030], [y + 0.012, y + 0.012], color=GRID,
            linewidth=1.2, zorder=2)
    y -= 0.014
    name, tag, makes = SCIENCE
    ax.text(LX, y, name, fontsize=11.5, color=INK, fontweight="bold",
            va="center", zorder=4)
    ax.text(TX, y, tag, fontsize=10.5, color=INK, va="center",
            fontfamily="monospace", zorder=4)
    ax.text(RX - 0.025, y, "→", fontsize=11.5, color=STEP, va="center",
            ha="center", zorder=4)
    ax.text(RX, y, makes, fontsize=11.5, color=INK, va="center", zorder=4)

    _finish(fig, ax, titles, path)


def _finish(fig, ax, titles, path):
    if titles:
        ax.text(0.0, 0.992, "Measuring KCRM red-channel throughput",
                fontsize=21, color=INK, fontweight="bold", va="top")
        ax.text(0.0, 0.938,
                "a standard star of known brightness, through the instrument, "
                "end to end",
                fontsize=13, color=INK, va="top")
    fig.subplots_adjust(left=0.026, right=0.988, top=0.988, bottom=0.012)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=ROOT)
    ap.add_argument("--pdf", action="store_true", help="also write the vector "
                    "copy for the deck")
    a = ap.parse_args()
    out = []
    # The talk copy: no calibration panel, no heading of its own, drawn as a
    # strip so it drops into a slide that already has a title.
    for name, kw in (
            ("kcrm_workflow", dict(detail=False, titles=False, figh=3.6)),
            ("kcrm_workflow_detail", dict(detail=True, titles=True, figh=7.5))):
        for ext in ("png",) + (("pdf",) if a.pdf else ()):
            path = os.path.join(a.outdir, f"{name}.{ext}")
            figure(path, **kw)
            out.append(os.path.basename(path))
    print(", ".join(out))


if __name__ == "__main__":
    main()
