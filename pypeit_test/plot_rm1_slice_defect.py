"""Detector column 653 is dead, and that is why 2024-04-01 setup B traced 25 slices.

RH1_PROCEDURE.md problem 4 found column 653 dying above spectral row ~3480 on
2023-11-08, and closed with "Column 653 may be dead on other KCRM data. If it is
a permanent detector [defect]...".  **It is.**  Measured here on 2024-04-01, a
different night, a different slicer and a different binning:

    median counts down each column, top eighth of the detector (rows 3612+)
        col 651   17788        col 653     290   <- dead
        col 652   17787        col 654   17550
                               col 655   17722

The interesting part is not that the column is dead.  It is that **the same dead
column breaks one setup of the night and not the other**, and the reason is
geometry: a dead column only invents a slice when it falls INSIDE a lit slice.
Where it falls in the gap between slices, it costs nothing, because a gap is
supposed to be dark there anyway.

    setup A (cenwave 6130)   at rows 3612+, the slice boundary already sits at
                             ~654, so column 653 lies on the edge and is
                             absorbed into it            -> 24 clean slices
    setup B (cenwave 7390)   at the same rows the slice is still fully lit out
                             to column 660, so column 653 sits ~8 columns INSIDE
                             it.  A drop from 17,800 to 290 is a far sharper
                             gradient than any real slit edge, so the tracer
                             registers an edge there     -> 25 slices, and
                             `run_pypeit` dies with "Alignment tracing has failed
                             on slit 4/25"

The grating setting moves the slice boundaries across the detector by a few
columns, so which of those two cases you get is decided by the central
wavelength.  That is why this cannot be predicted from the detector alone.

Fix: `[calibrations][[slitedges]] exclude_regions = 1:642:656,` -- trailing
comma required, because the parameter is a list.  Same mechanism as RH1's
`1:652:655,`, a slightly wider window.  Measured result: 24 slits, widths
139.9-140.7 px, rms 0.106, against an 8.9 px sliver and a 126.3 px truncated
neighbour before.  Note the window does clip the real boundary in the bottom
third of the detector, where the boundary leans through 642-656 and column 653
is still healthy; PypeIt traces the edge from the rows where it is visible and
carries it through, which was settled by running it rather than by argument.

Everything is drawn in **PypeIt's oriented detector frame**, which is the frame
the slit traces live in.  That is not a detail: `keck_kcrm` has `spatflip=True`,
so raw spatial column c becomes oriented column 4113-c, and the raw frame also
carries overscan that gets trimmed.  Plotting a raw frame under a slit trace
compares two different parts of the detector -- and hunting for the dead column
in raw coordinates finds nothing, because raw 653 is not oriented 653.

Vocabulary, used precisely throughout:

  spatial column   detector x, across the slices, AFTER trimming and the spatial
                   flip.  0-4113 on this Small/1x1 frame
  spectral row     detector y, along the dispersion.  0-4127
  slice            one of the 24 strips the image slicer cuts the field into.
                   ~140 spatial columns wide on Small/1x1
  slit edge        PypeIt's traced boundary of a slice: a left curve and a right
                   curve, each giving one spatial column per spectral row
  gap              the unilluminated strip between two slices
  sliver           a traced "slice" far narrower than 140 px.  Never physical

Usage:
    python plot_rm1_slice_defect.py [--snapshot slitsnap.npz] [-o out.png]
"""
import argparse
import glob
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.lines import Line2D                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "fits", "by_night", "RM1", "2024-04-01")
TRACE_B = "KR.20240402.16138.10.fits"      # setup B, cenwave 7390 -> 25 slits
TRACE_A = "KR.20240402.12636.17.fits"      # setup A, cenwave 6130 -> 24 slits

BLUE, ORANGE, AQUA, YELLOW, MAGENTA = ("#2a78d6", "#eb6834", "#1baf7a",
                                       "#eda100", "#e87ba4")
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"

C0, C1 = 610, 710
BANDS = [(0.03, 0.15), (0.28, 0.40), (0.53, 0.65), (0.78, 0.90), (0.93, 1.00)]
# Traced edges whose mean column lands in the window, measured from the
# pre-fix Slits: the first two have no cliff under them.
SPURIOUS = (646.66, 651.65)


def oriented(fname):
    """The flat as PypeIt sees it: overscan trimmed, spatial axis flipped."""
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.core import procimg
    spec = load_spectrograph("keck_kcrm")
    det, img, _, _, dsec, _ = spec.get_rawimage(os.path.join(RAW, fname), 1)
    return spec.orient_image(det, procimg.trim_frame(img, dsec < 1))


def edges_from(snapshot, key, setup):
    if snapshot and os.path.exists(snapshot):
        z = np.load(snapshot)
        return z[key + "_L"], z[key + "_R"]
    from pypeit.slittrace import SlitTraceSet
    g = glob.glob(os.path.join(ROOT, "RM1 pypeit run", "2024-04-01",
                               "pypeit_run", setup, "Calibrations", "Slits_*"))
    st = SlitTraceSet.from_file(g[0], chk_version=False)
    return np.asarray(st.left_init), np.asarray(st.right_init)


def show(ax, img, L, R, title, mark_spurious):
    sub = img[:, C0:C1]
    v1, v2 = np.percentile(sub[np.isfinite(sub)], [3, 99])
    ax.imshow(sub, origin="lower", aspect="auto", cmap="gray", vmin=v1, vmax=v2,
              extent=[C0, C1, 0, img.shape[0]], interpolation="nearest")
    rows = np.arange(L.shape[0])
    for i in range(L.shape[1]):
        for arr in (L, R):
            v = arr[:, i]
            m = np.nanmean(v)
            if not (C0 - 4 < m < C1 + 4):
                continue
            bad = mark_spurious and any(abs(m - s) < 2.0 for s in SPURIOUS)
            ax.plot(v, rows, color=ORANGE if bad else BLUE,
                    lw=2.6 if bad else 1.8, zorder=4 if bad else 3)
    ax.set_xlim(C0, C1)
    ax.set_ylim(0, L.shape[0])
    ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=9)
    ax.set_xlabel("spatial column (oriented detector x)", color=INK2, fontsize=9.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", default=None,
                    help="npz of the PRE-FIX slit traces.  Without it the "
                         "current Calibrations/ are read, which after the "
                         "exclude_regions fix no longer show the defect")
    ap.add_argument("-o", "--out",
                    default=os.path.join(ROOT, "rm1_7390_slice_defect.png"))
    args = ap.parse_args()

    ib, ia = oriented(TRACE_B), oriented(TRACE_A)
    Lb, Rb = edges_from(args.snapshot, "b7390", "keck_kcrm_B")
    La, Ra = edges_from(args.snapshot, "a6130", "keck_kcrm_A")

    fig = plt.figure(figsize=(15.2, 11.6), dpi=155)
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(2, 2, left=0.062, right=0.985, top=0.845, bottom=0.062,
                          hspace=0.30, wspace=0.20, height_ratios=[1.0, 0.88])
    axes = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (0, 1)]
    for ax in axes:
        ax.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTED)
            ax.spines[sp].set_linewidth(0.8)
        ax.tick_params(colors=INK2, labelsize=9, length=3)

    ny = ib.shape[0]

    # ---- 1  the symptom ----------------------------------------------------
    show(axes[0], ib, Lb, Rb,
         "1  setup B (cenwave 7390): the dome flat, with the edges PypeIt traced",
         True)
    axes[0].set_ylabel("spectral row (detector y)", color=INK2, fontsize=9.5)
    axes[0].axvline(653, color=AQUA, lw=1.2, ls="--", zorder=6)
    axes[0].annotate("column 653 — dead above row ~3600.\n"
                     "Here it sits ~8 columns inside a lit\n"
                     "slice, so the tracer reads it as an edge.",
                     xy=(653, 3900), xytext=(0.03, 0.955),
                     textcoords="axes fraction", color=AQUA, fontsize=9,
                     va="top", ha="left",
                     bbox=dict(boxstyle="round,pad=0.35", fc=SURFACE, ec="none",
                               alpha=0.9),
                     arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.5))
    axes[0].annotate("the real gap between two slices",
                     xy=(670, 1500), xytext=(0.40, 0.05),
                     textcoords="axes fraction", color=BLUE, fontsize=9,
                     va="bottom", ha="left",
                     bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec="none",
                               alpha=0.9),
                     arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.5))

    # ---- 2  the control ----------------------------------------------------
    show(axes[1], ia, La, Ra,
         "2  setup A (cenwave 6130), same night, same detector: 24 clean slices",
         False)
    axes[1].axvline(653, color=AQUA, lw=1.2, ls="--", zorder=6)
    axes[1].annotate("Column 653 is just as dead here (297 counts).\n"
                     "But at these rows the slice boundary already\n"
                     "sits at ~654, so the dead column lies ON the\n"
                     "edge and is absorbed into it. Same defect,\n"
                     "no spurious slice — the difference is where\n"
                     "the grating setting put the boundary.",
                     xy=(0.03, 0.955), xycoords="axes fraction", color=INK2,
                     fontsize=9, va="top", ha="left",
                     bbox=dict(boxstyle="round,pad=0.35", fc=SURFACE, ec="none",
                               alpha=0.9))

    # ---- 3  the cause, unambiguous -----------------------------------------
    ax = axes[2]
    eighth = ny // 8
    rows = np.array([(i + 0.5) * eighth for i in range(8)])
    for c, col, style in ((651, MUTED, "-"), (652, MUTED, "-"),
                          (654, MUTED, "-"), (655, MUTED, "-")):
        q = [np.median(ib[i * eighth:(i + 1) * eighth, c]) for i in range(8)]
        ax.plot(np.array(q) / 1000.0, rows, style, color=col, lw=1.5, alpha=0.75,
                zorder=3)
    q653 = [np.median(ib[i * eighth:(i + 1) * eighth, 653]) for i in range(8)]
    ax.plot(np.array(q653) / 1000.0, rows, "-o", color=AQUA, lw=2.6, ms=7,
            zorder=5, label="column 653")
    ax.plot([], [], color=MUTED, lw=1.5, label="columns 651, 652, 654, 655")
    ax.set_title("3  the cause: column 653 dies in the top eighth of the detector\n"
                 "median counts down each column, setup B dome flat",
                 color=INK, fontsize=10.5, loc="left", pad=9)
    ax.set_xlabel("median counts (thousands of ADU)", color=INK2, fontsize=9.5)
    ax.set_ylabel("spectral row (detector y)", color=INK2, fontsize=9.5)
    ax.grid(color=MUTED, alpha=0.18, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_ylim(0, ny)
    ax.annotate("290 ADU, against ~17,800 in\nevery neighbouring column",
                xy=(q653[-1] / 1000.0, rows[-1]), xytext=(6.0, 3050),
                color=AQUA, fontsize=9, va="center", ha="left",
                arrowprops=dict(arrowstyle="->", color=AQUA, lw=1.5))
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)

    # ---- 4  why one setup escapes ------------------------------------------
    ax = axes[3]
    r = np.arange(ny)
    dead_from = 3600
    ax.axhspan(dead_from, ny, color=AQUA, alpha=0.10, lw=0, zorder=0)
    ax.axvline(653, color=AQUA, lw=2.0, ls="--", zorder=4)
    ax.plot(Rb[:, 4], r, color=ORANGE, lw=2.2,
            label="setup B — slice boundary")
    ax.plot(Ra[:, 3], r, color=BLUE, lw=2.2,
            label="setup A — slice boundary")
    ax.set_title("4  why the same dead column breaks only setup B\n"
                 "where each setup's slice boundary sits, against column 653",
                 color=INK, fontsize=10.5, loc="left", pad=9)
    ax.set_xlabel("spatial column (oriented detector x)", color=INK2, fontsize=9.5)
    ax.set_ylabel("spectral row (detector y)", color=INK2, fontsize=9.5)
    ax.grid(color=MUTED, alpha=0.18, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_xlim(640, 672)
    ax.set_ylim(0, ny)
    ax.annotate("rows where column 653 is dead",
                xy=(641, (dead_from + ny) / 2), color=AQUA, fontsize=9,
                va="center", ha="left")
    ax.annotate("setup B: boundary is at ~662 here,\nso 653 is ~9 columns INSIDE\n"
                "the lit slice → spurious edge",
                xy=(662, 3850), xytext=(643.4, 2560), color=ORANGE, fontsize=9,
                va="top", ha="left",
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.5))
    ax.annotate("setup A: boundary is at ~655 here,\nso 653 is ON the edge → absorbed",
                xy=(655.6, 3850), xytext=(643.4, 1080), color=BLUE, fontsize=9,
                va="top", ha="left",
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.5))
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)

    fig.legend(handles=[Line2D([], [], color=BLUE, lw=2.2,
                               label="traced edge sitting on the gap (real)"),
                        Line2D([], [], color=ORANGE, lw=2.6,
                               label="traced edge with no gap under it (spurious)"),
                        Line2D([], [], color=AQUA, lw=2.0, ls="--",
                               label="detector column 653")],
               frameon=False, fontsize=10, loc="upper right",
               bbox_to_anchor=(0.985, 0.995), labelcolor=INK2, ncol=1)
    fig.suptitle("Detector column 653 is dead — and it only invents a slice "
                 "when it lands inside one",
                 color=INK, fontsize=15.5, x=0.006, ha="left", y=0.988)
    fig.text(0.006, 0.951,
             "2024-04-01, KCRM RM1, Small slicer 1x1 — why setup B traced 25 "
             "slices where KCRM has 24, and setup A the same night did not",
             color=INK2, fontsize=11, ha="left", va="top")
    fig.text(0.006, 0.916,
             "RH1_PROCEDURE.md problem 4 found column 653 dying on 2023-11-08 and "
             "asked whether it was permanent. It is: a different night, slicer and "
             "binning, same column, 290 ADU against ~17,800.\nWhat is new here is "
             "that the SAME dead column breaks one setup of a night and not the "
             "other. A dead column only invents a slice when it falls inside a lit "
             "slice; in a gap it costs nothing. The grating\nsetting moves the slice "
             "boundaries a few columns, and that decides which case you get — so it "
             "cannot be predicted from the detector alone. Fix: exclude_regions = "
             "1:642:656, — same mechanism as RH1's 1:652:655.",
             color=MUTED, fontsize=9.5, ha="left", va="top")
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
