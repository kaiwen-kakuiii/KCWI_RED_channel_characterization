#!/usr/bin/env python
"""Show which slices disagree with their own setup, and which slit is malformed.

Two failure modes, one figure:

LEFT  wavelength.  Every slice of a setup should start within ~40 A of its
      neighbours -- the slicer feeds the grating at slightly different angles.
      A slice whose template window was cut in the wrong place keeps a normal
      width and a smooth dispersion, so only this comparison exposes it.

RIGHT geometry.  Every KCRM slice is ~140 px wide.  A slice that is not is a
      tracing failure, and on 2023-11-08 the truncated one makes the astrometric
      transform return NaN, which kills the coadd.

Status colour marks state, never identity: grey is every healthy slice, red is a
slice that fails.  Both are labelled, so the colour never carries the meaning by
itself.

Usage:
    python plot_slice_diagnostics.py [-o out.png]
"""
import argparse
import glob
import os

import numpy as np
from pypeit.slittrace import SlitTraceSet
from pypeit.wavecalib import WaveCalib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RH1 pypeit run")

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8880"
CRITICAL, GRID = "#d03b3b", "#e6e5e0"
NORMAL = "#9a9a94"

# Measured before the 6520 setups were re-reduced.  Their WaveCalib files were
# deleted to force recalibration, so these four points cannot be re-derived from
# what is on disk -- they are quoted from the run that found them.
PRE_FIX = [("2023-10-17 C", 1310, -330.9), ("2023-10-17 C", 1458, -337.1),
           ("2023-10-20 B", 812, -336.6), ("2023-10-20 B", 1312, -201.2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out",
                    default=os.path.join(ROOT, "rh1_slice_diagnostics.png"))
    args = ap.parse_args()

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(13.2, 5.6), dpi=170, gridspec_kw=dict(width_ratios=[1.7, 1]))
    fig.patch.set_facecolor(SURFACE)

    # ---- LEFT: wavelength agreement, current state -------------------------
    axL.set_facecolor(SURFACE)
    setups, allofs = [], []
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        p = os.path.relpath(wc, RUNS).split(os.sep)
        tag = f"{p[0]} {p[2].replace('keck_kcrm_', '')}"
        wv = WaveCalib.from_file(wc, chk_version=False)
        blue = [np.asarray(wf.wave_soln).flatten().min()
                for wf in wv.wv_fits if wf is not None and wf.wave_soln is not None]
        if not blue:
            continue
        med = float(np.median(blue))
        offs = np.array(blue) - med
        setups.append(tag)
        allofs.append(offs)
        axL.plot(np.arange(len(offs)), offs, "o", ms=5, color=NORMAL,
                 alpha=0.65, zorder=3)

    # Three of the four sit within 6 A of each other, so their labels would
    # overlap; they are laddered instead and joined to their point by a leader.
    order = sorted(range(len(PRE_FIX)), key=lambda i: PRE_FIX[i][2])
    ladder = {}
    for rank, i in enumerate(order):
        ladder[i] = -352 + rank * 30
    for i, (tag, spat, off) in enumerate(PRE_FIX):
        axL.plot([24.8], [off], "X", ms=11, color=CRITICAL, zorder=5)
        axL.annotate(f"{tag}  spat {spat}   {off:+.0f} A",
                     xy=(24.8, off), xytext=(26.4, ladder[i]),
                     color=CRITICAL, fontsize=8.5, va="center",
                     arrowprops=dict(arrowstyle="-", color=CRITICAL,
                                     lw=0.8, alpha=0.6,
                                     shrinkA=0, shrinkB=4))

    axL.axhspan(-40, 40, color=NORMAL, alpha=0.14, zorder=1)
    axL.annotate("normal slice-to-slice spread (~40 A)", (0.4, 44),
                 color=MUTED, fontsize=8.5)
    axL.axhline(0, color=MUTED, lw=1, zorder=2)
    axL.set_xlabel("slice index within setup", color=INK2, fontsize=10)
    axL.set_ylabel("blue end minus setup median (A)", color=INK2, fontsize=10)
    axL.set_title("Wavelength: which slice disagrees with its setup",
                  color=INK, fontsize=12, loc="left", pad=10)
    axL.set_xlim(-1, 40)
    axL.set_ylim(-380, 70)
    axL.legend(handles=[
        Line2D([], [], ls="", marker="o", ms=6, color=NORMAL,
               label=f"all {sum(len(o) for o in allofs)} slices now on disk"),
        Line2D([], [], ls="", marker="X", ms=9, color=CRITICAL,
               label="shifted, measured before re-reduction")],
        frameon=False, fontsize=9, loc="lower left", labelcolor=INK2)

    # ---- RIGHT: slit geometry on the night that still fails -----------------
    axR.set_facecolor(SURFACE)
    sf = glob.glob(os.path.join(RUNS, "2023-11-08", "pypeit_run", "keck_kcrm_A",
                                "Calibrations", "Slits_*.fits.gz"))
    if sf:
        s = SlitTraceSet.from_file(sf[0], chk_version=False)
        L, R, _ = s.select_edges()
        w = np.median(R - L, axis=0)
        med = float(np.median(w))
        bad = np.abs(w - med) > 0.05 * med
        axR.bar(np.arange(len(w))[~bad], w[~bad], color=NORMAL, width=0.72, zorder=3)
        axR.bar(np.arange(len(w))[bad], w[bad], color=CRITICAL, width=0.72, zorder=4)
        for i in np.where(bad)[0]:
            axR.annotate(f"slice {i}\n{w[i]:.1f} px\nastrometric transform = NaN",
                         (i, w[i]), textcoords="offset points", xytext=(6, -34),
                         color=CRITICAL, fontsize=8.5)
        axR.axhline(med, color=MUTED, lw=1, ls="--", zorder=5)
        axR.annotate(f"median {med:.1f} px", (len(w) - 0.5, med),
                     textcoords="offset points", xytext=(-4, 5),
                     color=MUTED, fontsize=8.5, ha="right")
        axR.set_ylim(0, max(w) * 1.18)
    axR.set_xlabel("slice index", color=INK2, fontsize=10)
    axR.set_ylabel("slice width (detector px)", color=INK2, fontsize=10)
    axR.set_title("Geometry: 2023-11-08, the night that still fails",
                  color=INK, fontsize=12, loc="left", pad=10)

    for ax in (axL, axR):
        ax.grid(True, color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color("#d8d7d2")
        ax.tick_params(colors=INK2, labelsize=9)

    fig.tight_layout()
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"wrote {args.out}")
    print(f"  left: {len(setups)} setups, "
          f"{sum(len(o) for o in allofs)} slices, "
          f"max |offset| now {max(np.abs(np.concatenate(allofs))):.1f} A")


if __name__ == "__main__":
    main()
