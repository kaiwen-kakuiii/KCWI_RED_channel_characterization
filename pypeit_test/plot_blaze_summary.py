#!/usr/bin/env python
"""Combine every grating's throughput curves by BLAZE, and draw the superblaze.

A "blaze" here is one grating setting: a (grating, central wavelength, slicer)
triple.  Rotating a blazed grating to a new central wavelength moves its peak
efficiency with it, so two central wavelengths of one grating are two different
efficiency curves and must never be averaged together.  The slicer is part of
the key for a separate reason: a narrower slicer measures ~6-8% low (measured on
RL 2024-12-04, the only night carrying two deckers), so mixing deckers inside
one blaze would fold an instrumental offset into the scatter.

Everything else -- night, visit, standard star -- is a repeat of the same
measurement and IS combined.  Each cube is drawn thin over its own fitted range;
the mean is drawn thick over the range every cube in the group covers, because a
mean whose membership changes with wavelength puts steps in the curve that look
like structure.

Three outputs:

  blaze_plots/<grating>_blazes.png   one per grating, all its blazes
  kcrm_superblaze.png                per grating, the MAX over its blazes
  kcrm_composite.png                 per grating, the MEAN over its blazes
  blaze_summary.csv                  the numbers behind all of them
  trim_summary.csv                   how much wavelength each grating gave up
  edge_check.csv                     whether the 150 A edge rule is earned

Curves are read across the whole of SENS_ZEROPOINT_FIT_GPM.  Outside the GPM the
zeropoint polynomial extrapolates and has been seen to exceed 100%, so that stays
excluded; inside it, nothing is trimmed -- see EDGE_PROBE for why the old 150 A
edge rule was dropped.

Usage:
    python plot_blaze_summary.py [--outdir <dir>]
"""
import argparse
import collections
import csv
import glob
import os
import re

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import LSQUnivariateSpline

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402
import matplotlib.patheffects as pe                   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# The polynomial order each grating's results were settled at.  These are not
# interchangeable: order was tested per grating on a night pair and RL/RH1/RH2/
# RH4 landed on 15, RH3/RM1/RM2 on 5.  Reading the wrong suffix would silently
# plot a different reduction.
SUFFIX = {"RL": "IR_p15cov", "RH1": "IR_p15cov", "RH2": "IR_p15cov",
          "RH3": "IR_p5cov", "RH4": "IR_p15cov",
          "RM1": "IR_p5cov", "RM2": "IR_p5cov"}
ORDER = ["RL", "RM1", "RM2", "RH1", "RH2", "RH3", "RH4"]

# Cubes excluded from every mean, peak and envelope.  Both are RH1 nights whose
# photons did not arrive: compared within one star over 6050-6350 A with
# extinction already removed, they run 15% and 30% low in raw counts, so
# sensfunc and the zeropoint fit cannot be the cause.  Kept in the per-blaze
# figures in grey, so the figure shows what was removed.
DROP = {
    "g191b2b_2023-11-07_B": "non-photometric: 15% low in counts vs 2023-09-16, same star",
    "feige34_2024-03-13_B": "non-photometric: 30% low in counts vs 2023-11-17, same star",
    # RM2's only 9850 config, and its only coverage past ~9750 A.  Excluded on
    # the same standard as the two RH1 nights: its five frames fade 41% in nine
    # minutes while airmass AND seeing improve, and it reads 30.7% at 9600 A
    # where seven other configs read 36.2-40.1%.  Dropping it ends RM2's
    # envelope at 9748 A rather than carrying a cloud measurement into the red.
    "feige110_2023-09-23_B": "non-photometric: fades 41% in 9 min while airmass and seeing improve",
}

# Kept, but a lower limit -- do not read these as the grating's efficiency.
FLAG = {
    "feige34_2024-03-15_B": "-6.9% against its own grating for reasons this run does not resolve",
}

# Validated categorical palette, 7 slots, assigned in fixed order -- to the
# grating in the superblaze, to the night in a per-blaze figure.  Checks on the
# light surface #fcfcfb: worst adjacent CVD dE 9.1 (protan), worst adjacent
# normal-vision dE 19.6, contrast WARNs on aqua/yellow/pink.  That WARN
# obligates relief, which is why every figure here carries direct labels and
# every number also appears in GRATING_SUMMARY.md.
# Matched to the KCRM predicted-throughput figure this run is compared against,
# slot for slot in ORDER, so the two can be read side by side without a mental
# remap.  This is a deliberate override of the validated palette: it FAILS CVD
# separation on RH1 green vs RH2 orange (dE 5.7 protan, floor 6) and drops below
# 3:1 contrast on RM2 yellow (1.8), RH2 orange (2.75) and RH3 teal (2.19).
# Comparability against the reference figure was judged worth it; the relief is
# that every curve on the all-grating figures is directly labelled AND legended,
# and every number also appears in GRATING_SUMMARY.md.  RH1 and RH2 overlap in
# wavelength, so that pair is the one to watch on a projector.
HUES = ["#4285F4", "#DB4437", "#F4B400", "#0F9D58", "#FF6D00", "#46BDC6",
        "#AB30C4"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SPINE, SURFACE = "#e6e5e0", "#d8d7d2", "#fcfcfb"
# The 150 A fit-edge rule this analysis used to apply is GONE, removed
# 2026-08-27.  It was inherited from the RL run, where it was justified by
# comparing an interior RMS (0.0300 mag) against a worst-single-point edge error
# (0.117 mag) -- two different statistics, so their ratio never measured edge
# degradation.  Recomputing both statistics in both zones on the current
# coverage-trimmed products put the edge/interior RMS ratio at 0.94-1.39 and the
# max ratio at 0.93-1.16: the zeropoint fit is as good at its boundary as in its
# middle.  The zone the rule protected against -- partial slice coverage -- is
# now removed BEFORE the fit by `measure_trim --spec2d`, so the rule was doing
# that job a second time, at a cost of a flat 300 A per config (47% of RH1's
# fitted range).  Everything below now reads the WHOLE fitted range.
#
# `edge_residuals` still probes at 150 A; that is a diagnostic, not a cut, and it
# is what keeps the decision falsifiable.
EDGE_PROBE = 150.0

# Knot spacing for the spline that smooths both all-grating figures.  600 A was
# chosen by fitting 400/600/900 against the unsmoothed envelopes: 400 was still
# following RH1's hand-over dip, 900 was stiff enough to miss RM1's blue rise
# and flatten RM2's 8000 A shoulder.  One value serves all seven gratings --
# unlike a window, a knot spacing does not have to be retuned per span.
KNOT_A = 600.0
MIN_KNOTS = 2

# Where each grating's direct label hangs off its own envelope.  Seven envelopes
# overlap heavily between 6000 and 7500 A, so "at the maximum, offset upward"
# collides; these were placed by looking at the render.  (wavelength anchor or
# None for the envelope max, dx, dy in points, ha)
LABEL = {"RL":  (7600, 0, -20, "center"),
         "RM1": (None, 0, 16, "center"),
         "RM2": (None, 0, 16, "center"),
         "RH1": (6200, -34, 12, "right"),
         "RH2": (6980, 0, 15, "center"),
         "RH3": (None, 0, 16, "center"),
         "RH4": (None, 30, 12, "left")}


def read_cube(sens, grating):
    """One reduced cube: its interior throughput curve and its configuration."""
    sdir = os.path.dirname(sens)
    tag = os.path.basename(sens).split("_sens_")[0]
    with fits.open(sens) as h:
        t = Table(h["SENS"].data)
        wave = np.asarray(h["WAVE"].data).flatten()
        thru = np.asarray(h["THROUGHPUT"].data).flatten() * 100.0
    gpm = np.asarray(t["SENS_ZEROPOINT_FIT_GPM"][0], dtype=bool)
    sw = np.asarray(t["SENS_WAVE"][0]).flatten()
    lo, hi = float(sw[gpm].min()), float(sw[gpm].max())
    # SENS_WAVE is zero-padded to a fixed length; its non-zero part is what
    # pypeit_sensfunc was actually handed, i.e. the spectrum AFTER
    # measure_trim.py's counts and slice-coverage trim.
    nz = sw > 1
    sens_lo, sens_hi = float(sw[nz].min()), float(sw[nz].max())
    ok = (wave >= lo) & (wave <= hi) & np.isfinite(thru) & (thru > 0)
    if not ok.any():
        return None
    w, v = wave[ok], thru[ok]

    meta = {}
    pf = glob.glob(os.path.join(sdir, "*.pypeit"))
    sci = []
    if pf:
        txt = open(pf[0]).read()
        for k in ("cenwave", "decker", "binning"):
            m = re.search(rf"^\s*{k}:\s*(\S+)", txt, re.M)
            if m:
                meta[k] = m.group(1)
        for line in txt.splitlines():
            if ".fits" in line and "|" in line:
                f = [c.strip() for c in line.split("|")]
                if len(f) > 2 and f[1] == "science":
                    sci.append(f)
    # Round the central wavelength: the same nominal setting is written
    # 6519.939 one night and 6519.996 another, which would split one blaze in two.
    cw = round(float(meta["cenwave"])) if "cenwave" in meta else None

    # Frames actually coadded, not frames available: the drivers drop sub-10 s
    # acquisition frames and (RM2) starless ones, and the .coadd3d is the record
    # of what the cube is really made of.
    used = []
    # Name the coadd3d after the cube where one exists.  A setup directory can
    # hold several -- RH1's two airmass visits, RH4's two standards -- and
    # taking the first match would report one cube's frame list against the
    # other cube's curve.  RL and RH2 wrote a single <star>.coadd3d instead, so
    # fall back to a lone file, never to an ambiguous one.
    co = [os.path.join(sdir, f"{tag}.coadd3d")]
    if not os.path.exists(co[0]):
        co = glob.glob(os.path.join(sdir, "*.coadd3d"))
        if len(co) != 1:
            co = []
    if co:
        blk = open(co[0]).read().split("spec2d read")[-1].split("spec2d end")[0]
        used = [l.strip() for l in blk.splitlines()
                if l.strip() and l.strip() != "filename"]
    exps = []
    for u in used:
        m = re.search(r"spec2d_(KR\.[0-9.]+)-", os.path.basename(u))
        if not m:
            continue
        for f in sci:
            if m.group(1) in f[0]:
                try:
                    exps.append(float(f[-8]))
                except (ValueError, IndexError):
                    pass
                break

    # What the reduction delivered, before any trim: the extracted spectrum's
    # own extent, which equals the datacube's wavelength axis (extraction cuts
    # spatially, never in wavelength).
    s1 = os.path.join(sdir, f"spec1d_{tag}.fits")
    deliv_lo = deliv_hi = np.nan
    if os.path.exists(s1):
        w1 = np.asarray(fits.getdata(s1, 1)["BOX_WAVE"]).flatten()
        w1 = w1[w1 > 1]
        if w1.size:
            deliv_lo, deliv_hi = float(w1.min()), float(w1.max())

    cube = os.path.join(sdir, f"{tag}.fits")
    am = float(fits.getheader(cube).get("AIRMASS", np.nan)) \
        if os.path.exists(cube) else np.nan

    # Every point inside the fit is used.  Kept as a mask rather than deleted so
    # the shape of the code still says "this is the usable part of the curve",
    # which is where a future cut would go if one is ever earned.
    interior = np.ones(w.shape, dtype=bool)
    peak = pw = mean = np.nan
    if interior.any():
        i = int(np.argmax(v[interior]))
        peak, pw = float(v[interior][i]), float(w[interior][i])
        mean = float(v[interior].mean())
    night = os.path.basename(os.path.dirname(os.path.dirname(sdir)))
    return dict(grating=grating, tag=tag, night=night,
                setup=os.path.basename(sdir).replace("keck_kcrm_", ""),
                star=tag.split("_")[0].lower(), cw=cw,
                decker=meta.get("decker", "?"), binning=meta.get("binning", "?"),
                nframes=len(used), exptime=sum(exps), airmass=am,
                lo=lo, hi=hi, wave=w, thru=v, interior=interior,
                deliv_lo=deliv_lo, deliv_hi=deliv_hi,
                sens_lo=sens_lo, sens_hi=sens_hi,
                peak=peak, peak_wave=pw, mean=mean,
                drop=DROP.get(tag), flag=FLAG.get(tag))


def collect():
    out = []
    for g in ORDER:
        pat = os.path.join(ROOT, f"{g} pypeit run", "*", "pypeit_run",
                           "keck_kcrm_*", f"*_sens_{SUFFIX[g]}.fits")
        for sens in sorted(glob.glob(pat)):
            if "order_study" in sens:
                continue
            c = read_cube(sens, g)
            if c:
                out.append(c)
    return out


def combine(group, step=1.0):
    """Mean curve over the range EVERY member covers, on a 1 A grid.

    Restricting to the common range is deliberate.  Averaging wherever at least
    one cube contributes makes the membership change with wavelength, and each
    change puts a step in the mean that reads as instrument structure.
    """
    keep = [c for c in group if not c["drop"]]
    if not keep:
        return None, None, None
    lo = max(c["lo"] for c in keep)
    hi = min(c["hi"] for c in keep)
    if hi - lo < 50:                    # no usable common range
        return None, None, None
    grid = np.arange(lo, hi + step, step)
    stack = np.vstack([np.interp(grid, c["wave"], c["thru"]) for c in keep])
    return grid, stack.mean(axis=0), stack.std(axis=0)


def smooth(x, y, knot_A=KNOT_A):
    """Least-squares cubic spline with knots every `knot_A` Angstroms.

    Replaces the Savitzky-Golay filter that was here before, because a moving
    window is the wrong tool for this curve.  What has to go is the SCALLOP --
    the dip where the envelope hands over from one central-wavelength setting to
    the next, or where a blaze enters the composite average.  That is an
    artifact of how the grating was scheduled, not of the optics, and its width
    is set by the gap between settings, which differs per grating.

    A window wide enough to erase RH1's 320 A hand-over is a large fraction of
    RH1's whole 1090 A range, and a filter whose window approaches the length of
    its data returns a polynomial, not a measurement.  A spline has no such
    failure mode: stiffness comes from knot SPACING, so the same 600 A setting
    is gentle on RL's 4273 A envelope and firm on RH1's 1090 A one, and both
    stay faithful.

    Measured against the unsmoothed curves: this tracks RL's real structure --
    the 9200 A peak and the 7000-8500 A plateau -- more closely than the filter
    did, while removing RH1's double bump, RH2's shoulder and RH3's kink
    entirely.  It leaves RM2 and RH4 essentially untouched, which is the right
    answer for the two gratings whose settings are too closely spaced (RM2) or
    too few (RH4, one blaze) to scallop in the first place.

    A floor of MIN_KNOTS interior knots applies, because the five narrower
    gratings span under 1200 A and a 600 A spacing would give them none at all
    -- a bare quartic, which cannot follow an asymmetric blaze.  Measured on the
    envelopes, dropping to zero knots costs RH1 1.03 points RMS against its
    unsmoothed curve; two knots brings that to a level comparable with the other
    gratings while still returning a single smooth hump.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 8:
        return y
    nknot = max(MIN_KNOTS, int((x[-1] - x[0]) // knot_A) - 1)
    nknot = min(nknot, max(0, x.size // 4 - 2))
    if nknot <= 0:
        return np.polyval(np.polyfit(x, y, min(4, x.size - 1)), x)
    t = np.linspace(x[0], x[-1], nknot + 2)[1:-1]
    return LSQUnivariateSpline(x, y, t, k=3)(x)


def peak_stats(w, v, tol=0.5):
    """Peak, its wavelength, and how wide the top is.

    The plateau -- the span over which the curve stays within `tol` points of
    its maximum -- is reported because argmax alone is misleading on a broad
    grating.  RL's curve is flat to a few tenths of a point over thousands of
    Angstroms and carries three near-equal local maxima, so its "peak
    wavelength" moves between them on noise; RH4's top is a few hundred A wide
    and means something.  The number tells the reader which case they have.
    """
    i = int(np.argmax(v))
    within = np.flatnonzero(v >= v[i] - tol)
    return float(v[i]), float(w[i]), float(w[within[0]]), float(w[within[-1]])


def ylow(values, pct=5.0, pad=1.5):
    """Bottom of the y-axis: low enough to show the curves, not the tails.

    With the 150 A edge rule gone, every curve now runs all the way to its fit
    boundary, and a fit boundary is where throughput is falling off a cliff --
    RH4 reaches 13.7%, RM1 8.7%.  Scaling to those terminal tails squeezes the
    25-41% band the figures exist to compare into the top half of the axes.
    Cutting at the 5th percentile of the plotted values clips ~5% of points,
    all of it inside the last few tens of Angstroms of some config, and the
    captions say the tails run off the bottom so a clipped line is not read as
    data ending.
    """
    v = np.asarray(values)
    v = v[np.isfinite(v)]
    return max(0.0, float(np.percentile(v, pct)) - pad)


def style_axes(ax, xlabel="wavelength (A)", ylabel="throughput (%)",
               label_fs=10, tick_fs=9):
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(SPINE)
    ax.tick_params(colors=INK2, labelsize=tick_fs)
    ax.set_xlabel(xlabel, color=INK2, fontsize=label_fs)
    ax.set_ylabel(ylabel, color=INK2, fontsize=label_fs)


def grating_figure(g, gg, path):
    """Every blaze of one grating on one axis, nights combined within a blaze.

    Colour keys on the BLAZE -- the (cenwave, slicer) setting -- because that is
    the thing that makes two curves comparable or not.  Night, airmass visit and
    standard star are combined into the blaze mean rather than given hues of
    their own; they survive as the faint traces behind each mean and as the
    +/-1 sd band around it, which is what night-to-night agreement looks like.

    `gg` is {(grating, cenwave, decker): [cubes]} restricted to this grating.
    """
    keys = sorted(gg, key=lambda k: (k[1], k[2]))
    colour = {k: HUES[i % len(HUES)] for i, k in enumerate(keys)}
    assert len(keys) <= len(HUES), f"{g}: {len(keys)} blazes, {len(HUES)} hues"

    fig, ax = plt.subplots(figsize=(11.4, 6.0), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    every, drawn, any_drop = [], [], False
    for k in keys:
        group = gg[k]
        every.extend(group)
        kept = [c for c in group if not c["drop"]]
        col = colour[k]
        any_drop = any_drop or len(kept) < len(group)

        # Individual cubes sit behind the mean as texture, not as series: with
        # eight cubes in RH1's 6520 alone, drawing them at series weight would
        # bury the six other blazes this figure exists to compare.
        for c in group:
            cc = MUTED if c["drop"] else col
            m = c["interior"]
            if m.any():
                ax.plot(c["wave"][m], c["thru"][m], color=cc, lw=0.9,
                        alpha=0.15 if not c["drop"] else 0.5, zorder=3)
        if not kept:
            continue

        grid, mean, sd = combine(group)
        if grid is None:
            c = kept[0]
            m = c["interior"]
            if not m.any():
                continue
            grid, mean, sd = c["wave"][m], c["thru"][m], None
        if sd is not None and len(kept) > 1:
            ax.fill_between(grid, mean - sd, mean + sd, color=col, alpha=0.16,
                            lw=0, zorder=4)
        ax.plot(grid, mean, color=col, lw=2.6, zorder=5,
                solid_capstyle="round")
        pk, pw, pl, ph = peak_stats(grid, mean)
        ax.plot([pw], [pk], "o", ms=7.5, mfc=SURFACE, mec=col, mew=2, zorder=7)
        drawn.append((k, col, len(kept), len({c["night"] for c in kept}),
                      pk, pw))

    # Direct labels only while they fit: past four they collide with the curves
    # they name, and the legend carries identity instead.
    if len(drawn) <= 4:
        for k, col, ncube, nnight, pk, pw in drawn:
            # Hung to the RIGHT of the marker, not above it: adjacent blazes
            # peak within tens of Angstroms of each other but a point or two
            # apart in throughput, so the vertical gap separates the labels
            # where the horizontal one does not.
            ax.annotate(f"{k[1]} A", (pw, pk), textcoords="offset points",
                        xytext=(11, -1), ha="left", va="center", color=col,
                        fontsize=10,
                        fontweight="bold", zorder=8,
                        path_effects=[pe.withStroke(linewidth=3.5,
                                                    foreground=SURFACE)])

    style_axes(ax)
    sub = ("one line per blaze -- the (central wavelength, slicer) setting"
           "   ·   thick = mean of that blaze's cubes, band = +/-1 sd,"
           " thin = the cubes themselves\n"
           "circle = blaze peak   ·   the whole fitted range is shown; the "
           "steepest fit-end tails run off the bottom")
    if any_drop:
        sub += "   ·   grey = excluded, non-photometric"
    nsub = sub.count("\n") + 1
    ax.set_title(f"KCRM {g} throughput, one curve per blaze",
                 color=INK, fontsize=13.5, loc="left", pad=18 + 15 * nsub)
    ax.annotate(sub, xy=(0, 1.0), xycoords="axes fraction",
                textcoords="offset points", xytext=(0, 16), color=MUTED,
                fontsize=9, va="bottom", ha="left", annotation_clip=False)

    handles = []
    for k, col, ncube, nnight, pk, pw in drawn:
        handles.append(Line2D([], [], color=col, lw=2.6,
                              label=f"{k[1]} A · {k[2]} · {nnight} night"
                                    f"{'s' if nnight != 1 else ''}, "
                                    f"{ncube} cube{'s' if ncube != 1 else ''}"
                                    f"  —  {pk:.1f}% @ {pw:.0f} A"))
    if any_drop:
        handles.append(Line2D([], [], color=MUTED, lw=1.6,
                              label="excluded night"))
    leg = ax.legend(handles=handles, frameon=False, fontsize=9,
                    loc="center left", bbox_to_anchor=(1.005, 0.5),
                    labelcolor=INK2)
    leg.set_title("blaze", prop={"size": 9})
    leg.get_title().set_color(MUTED)

    # y-range from the SURVIVING cubes only.  A non-photometric night can sit 10
    # points low, and letting it set the floor squeezes every real curve into
    # the top third of the axes; the grey trace clips instead, which is the
    # right priority for context that is already labelled as excluded.
    iv = np.concatenate([c["thru"][c["interior"]] for c in every
                         if c["interior"].any() and not c["drop"]])
    ax.set_ylim(ylow(iv), iv.max() + 3.0)
    ax.set_xlim(min(c["lo"] for c in every) - 30,
                max(c["hi"] for c in every) + 30)
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def union(intervals):
    """Merge intervals; return (merged list, total span)."""
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out, sum(b - a for a, b in out)


def trim_report(cubes, csv_path):
    """How much wavelength each grating gave up, and at which of the three cuts.

    A config's spectrum is narrowed twice, and only the first cut is large:

      1. `measure_trim.py`, before sensfunc -- drops the ends on a counts
         criterion AND on full slice coverage, because near the cube ends only
         some of the 24 slices contribute and the coadd fabricates flux there.
      2. the zeropoint fit itself -- SENS_ZEROPOINT_FIT_GPM is usually the whole
         input, a few Angstroms less where the fit rejected an end.

    A third cut, this analysis' 150 A fit-edge rule, was removed 2026-08-27; see
    EDGE_PROBE at the top of this file for the measurement that retired it.

    Both the per-config mean and the UNION over a grating's configs are
    reported.  The union is what the grating actually covers once every central
    wavelength is stacked.
    """
    rows = []
    for g in ORDER:
        cc = [c for c in cubes if c["grating"] == g and not c["drop"]
              and np.isfinite(c["deliv_lo"])]
        if not cc:
            continue
        deliv = np.array([c["deliv_hi"] - c["deliv_lo"] for c in cc])
        sens = np.array([c["sens_hi"] - c["sens_lo"] for c in cc])
        fit = np.array([c["hi"] - c["lo"] for c in cc])
        inter = fit
        du, ds = union([[c["deliv_lo"], c["deliv_hi"]] for c in cc])
        iu, isp = union([[c["lo"], c["hi"]] for c in cc])
        gaps = [(iu[i][1], iu[i + 1][0]) for i in range(len(iu) - 1)]
        rows.append(dict(
            grating=g, ncubes=len(cc), deliv=deliv.mean(), sens=sens.mean(),
            fit=fit.mean(), interior=inter.mean(),
            cut_trim=(deliv - sens).mean(), cut_fit=(sens - fit).mean(),
            kept=100 * inter.mean() / deliv.mean(),
            u_lo=du[0][0], u_hi=du[-1][1], u_span=ds,
            ui_lo=iu[0][0], ui_hi=iu[-1][1], ui_span=isp,
            u_kept=100 * isp / ds,
            gaps="; ".join(f"{a:.0f}-{b:.0f}" for a, b in gaps)))

    with open(csv_path, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "ncubes", "delivered_A", "after_trim_std_A",
                      "fit_range_A", "cut_trim_std_A", "cut_fit_A",
                      "kept_pct",
                      "union_lo_A", "union_hi_A", "union_span_A",
                      "union_fit_lo_A", "union_fit_hi_A",
                      "union_fit_span_A", "union_kept_pct",
                      "interior_gaps_A"])
        for r in rows:
            wtr.writerow([r["grating"], r["ncubes"], f"{r['deliv']:.0f}",
                          f"{r['sens']:.0f}", f"{r['fit']:.0f}",
                          f"{r['cut_trim']:.0f}", f"{r['cut_fit']:.0f}",
                          f"{r['kept']:.1f}", f"{r['u_lo']:.0f}",
                          f"{r['u_hi']:.0f}", f"{r['u_span']:.0f}",
                          f"{r['ui_lo']:.0f}", f"{r['ui_hi']:.0f}",
                          f"{r['ui_span']:.0f}", f"{r['u_kept']:.1f}",
                          r["gaps"]])
    return rows


def edge_residuals(csv_path):
    """Is the zeropoint fit actually worse near its boundary?  Measure it.

    This is the measurement that retired the 150 A rule, kept running so the
    decision stays falsifiable rather than becoming folklore.

    THROUGHPUT_PROCEDURE justified the rule with "~3% interior, up to ~12%
    within 150 A of a boundary".  Those two numbers are an RMS (interior
    scatter, 0.0300 mag) and a MAXIMUM (worst edge error, 0.117 mag) --
    different statistics, so the ratio between them never measured edge
    degradation.

    This recomputes both statistics in both zones, per grating, on the current
    coverage-trimmed products.  It reads |SENS_ZEROPOINT - SENS_ZEROPOINT_FIT|
    inside the fit, bins it by distance from the nearer boundary, and writes the
    edge (<EDGE_PROBE) and interior (>2*EDGE_PROBE) RMS and max.  If a future
    reduction ever makes the edge genuinely worse, the ratios here move first.

    Note the zones are absolute Angstroms, which is the rule as it was written.
    RH1's whole fitted range is 554 A, so it has NO interior zone by that
    definition -- every wavelength it measures sat inside what a rule derived on
    RL called an edge.  That alone should have flagged the rule.
    """
    rows = []
    for g in ORDER:
        pat = os.path.join(ROOT, f"{g} pypeit run", "*", "pypeit_run",
                           "keck_kcrm_*", f"*_sens_{SUFFIX[g]}.fits")
        dist, resid, spans = [], [], []
        for sens in sorted(glob.glob(pat)):
            if "order_study" in sens:
                continue
            tag = os.path.basename(sens).split("_sens_")[0]
            if tag in DROP:
                continue
            with fits.open(sens) as h:
                t = Table(h["SENS"].data)
            sw = np.asarray(t["SENS_WAVE"][0]).flatten()
            zp = np.asarray(t["SENS_ZEROPOINT"][0]).flatten()
            zf = np.asarray(t["SENS_ZEROPOINT_FIT"][0]).flatten()
            gpm = np.asarray(t["SENS_ZEROPOINT_FIT_GPM"][0], dtype=bool)
            m = gpm & np.isfinite(zp) & np.isfinite(zf) & (zp > 0)
            if m.sum() < 50:
                continue
            w, r = sw[m], np.abs(zp[m] - zf[m])
            lo, hi = w.min(), w.max()
            dist.append(np.minimum(w - lo, hi - w))
            resid.append(r)
            spans.append(hi - lo)
        if not resid:
            continue
        d, r = np.concatenate(dist), np.concatenate(resid)
        e, i = d < EDGE_PROBE, d >= 2 * EDGE_PROBE
        row = dict(grating=g, span=float(np.mean(spans)), npix=int(r.size),
                   e_rms=float(np.sqrt((r[e] ** 2).mean())) if e.sum() > 20 else np.nan,
                   e_max=float(r[e].max()) if e.sum() else np.nan,
                   i_rms=float(np.sqrt((r[i] ** 2).mean())) if i.sum() > 20 else np.nan,
                   i_max=float(r[i].max()) if i.sum() > 20 else np.nan)
        rows.append(row)
    with open(csv_path, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "mean_fit_span_A", "npoints",
                      "edge_rms_mag", "edge_max_mag",
                      "interior_rms_mag", "interior_max_mag",
                      "rms_ratio", "max_ratio"])
        for r in rows:
            rr = ("" if not np.isfinite(r["i_rms"]) else f"{r['e_rms']/r['i_rms']:.2f}")
            mr = ("" if not np.isfinite(r["i_max"]) else f"{r['e_max']/r['i_max']:.2f}")
            wtr.writerow([r["grating"], f"{r['span']:.0f}", r["npix"],
                          f"{r['e_rms']:.3f}", f"{r['e_max']:.3f}",
                          "" if not np.isfinite(r["i_rms"]) else f"{r['i_rms']:.3f}",
                          "" if not np.isfinite(r["i_max"]) else f"{r['i_max']:.3f}",
                          rr, mr])
    return rows


def blaze_means(gg):
    """[(slicer, wavelengths, mean, is_lower_limit)] for one grating's blazes.

    Same curves the envelope and the composite are both built from, so the two
    all-grating figures can never be drawn off different inputs.
    """
    curves = []
    for (gr, cw, deck), grp in sorted(gg.items(), key=lambda kv: kv[0][1]):
        kept = [c for c in grp if not c["drop"]]
        if not kept:
            continue
        grid, mean, _ = combine(grp)
        if grid is None:          # one cube, or no common range: use it alone
            c = kept[0]
            m = c["interior"]
            if not m.any():
                continue
            grid, mean = c["wave"][m], c["thru"][m]
        curves.append((deck, np.asarray(grid), np.asarray(mean),
                       all(c["flag"] for c in kept)))
    return curves


# Fraction of its own width over which a blaze fades into and out of the
# composite, and the hard cap in Angstroms.  A blaze that switched on at full
# weight put a STEP in the mean at the wavelength it entered -- and it entered
# on its own steep flank, at the edge of its own fitted range, which is the
# least trustworthy part of it.  Ramping the weight there removes the step at
# its source instead of smearing it out with a wider filter.  15% of the blaze's
# own width keeps the ramp proportionate (RH1's ~550 A blazes get ~80 A, RL's
# ~3300 A ones would get 500 and are capped at 250).
TAPER_FRAC, TAPER_MAX = 0.15, 250.0

def mean_over_blazes(curves, knot_A=KNOT_A):
    """Weighted mean of blaze curves on a 1 A grid, tapered at each blaze's ends.

    Returns (grid, mean, sd, n, neff): `n` counts blazes with any weight at all,
    `neff` is the summed weight -- the effective number contributing, which is
    what the mean is actually an average of.
    """
    lo = min(c[1].min() for c in curves)
    hi = max(c[1].max() for c in curves)
    grid = np.arange(lo, hi + 1.0, 1.0)
    stack, wts = [], []
    for _, w, v, _ in curves:
        y = np.interp(grid, w, v, left=np.nan, right=np.nan)
        t = min(TAPER_MAX, TAPER_FRAC * (w[-1] - w[0]))
        if t <= 0:
            wt = np.isfinite(y).astype(float)
        else:
            # Raised cosine: 0 at the blaze's own edge, 1 by `t` inside it.
            d = np.minimum(grid - w[0], w[-1] - grid)
            wt = 0.5 * (1 - np.cos(np.pi * np.clip(d / t, 0, 1)))
            wt[~np.isfinite(y)] = 0.0
        stack.append(y)
        wts.append(wt)
    stack, wts = np.vstack(stack), np.vstack(wts)
    n = (wts > 0).sum(axis=0)
    neff = wts.sum(axis=0)
    ok = neff > 0
    val = np.where(np.isfinite(stack), stack, 0.0)
    mean = np.full(grid.shape, np.nan)
    mean[ok] = (val[:, ok] * wts[:, ok]).sum(axis=0) / neff[ok]
    sd = np.zeros(grid.shape)
    multi = n > 1
    m2 = multi & ok
    if m2.any():
        d2 = (val[:, m2] - mean[m2]) ** 2
        sd[m2] = np.sqrt((d2 * wts[:, m2]).sum(axis=0) / neff[m2])
    mean[ok] = smooth(grid[ok], mean[ok], knot_A)
    if m2.any():
        sd[m2] = smooth(grid[m2], sd[m2], knot_A)
    return grid, mean, sd, n, neff


def slicer_drag(curves):
    """How far the Small-slicer blazes pull a composite down, in points.

    Worth reporting separately from the envelope's treatment of the same
    blazes.  In a MAX a Small blaze only matters where it wins outright; in a
    MEAN it is in every average it touches, so the ~6-8% narrow-slicer penalty
    leaks into the composite wherever a Small config overlaps.  Returns
    (mean, worst) over the range the non-Small blazes cover, or None when the
    grating has no Small config or nothing else to compare against.
    """
    small = [c for c in curves if c[0] == "Small"]
    rest = [c for c in curves if c[0] != "Small"]
    if not small or not rest:
        return None
    g_all, m_all, _, _, _ = mean_over_blazes(curves)
    g_big, m_big, _, _, _ = mean_over_blazes(rest)
    d = np.interp(g_big, g_all, m_all) - m_big
    return float(np.nanmean(d)), float(np.nanmin(d))


def composite(groups, path, csv_path):
    """Per grating, the MEAN of its blazes at each wavelength, not the max.

    The envelope answers "what can this grating reach here, given the right
    grating angle".  This answers a different question -- "what did this grating
    deliver here, averaged over the angles actually observed" -- and the two
    differ by more than presentation.

    Away from a blaze peak every contributing curve is on its own falling flank,
    so the composite is pulled BELOW what the grating can do at that wavelength,
    by an amount that depends only on which central wavelengths happened to be
    scheduled.  It is the honest summary of this dataset and the wrong number to
    quote as the grating's efficiency; the envelope is the right one.  Both are
    drawn from the identical blaze means (`blaze_means`), so the gap between the
    two figures is exactly that sampling effect and nothing else.

    Drawn as one solid line per grating, with no scatter band and no marking of
    the stretches a single blaze carries.  The +/-1 sd is still computed and
    still written to composite_by_grating.csv alongside the blaze count -- it is
    off the figure, not out of the data.
    """
    fig, ax = plt.subplots(figsize=(13.2, 6.6), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    colour = {g: HUES[i] for i, g in enumerate(ORDER)}
    rows, summary, drawn_y = [], [], []

    for g in ORDER:
        curves = blaze_means({k: v for k, v in groups.items() if k[0] == g})
        if not curves:
            continue
        grid, mean, sd, n, neff = mean_over_blazes(curves)
        ok, multi = neff > 0, n > 1
        # One solid line per grating, and no scatter band.  The stretches
        # carried by a single blaze are not marked either: the figure is for
        # comparing seven gratings at a glance, and the per-wavelength blaze
        # count is a column of composite_by_grating.csv for anyone who needs it.
        ax.plot(grid, mean, color=colour[g], lw=2.8, zorder=4,
                solid_capstyle="round")
        drawn_y.append(mean[ok])
        col = colour[g]

        anchor, dx, dy, ha = LABEL.get(g, (None, 0, 16, "center"))
        if anchor is None or not (grid[0] <= anchor <= grid[-1]):
            j = int(np.nanargmax(mean))
        else:
            j = int(np.argmin(np.abs(grid - anchor)))
        ax.annotate(g, (grid[j], mean[j]), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, color=col, fontsize=12.5,
                    fontweight="bold", zorder=8,
                    path_effects=[pe.withStroke(linewidth=3.5,
                                                foreground=SURFACE)])
        for i in np.flatnonzero(ok):
            rows.append((g, grid[i], mean[i], sd[i] if multi[i] else np.nan,
                         int(n[i]), neff[i]))
        i = int(np.nanargmax(mean))
        drag = slicer_drag(curves)
        summary.append(dict(grating=g, peak=mean[i], peak_wave=grid[i],
                            n_min=int(n[ok].min()), n_max=int(n.max()),
                            nblaze=len(curves),
                            drag_mean=drag[0] if drag else None,
                            drag_worst=drag[1] if drag else None))

    style_axes(ax)
    ax.set_title("KCRM red channel: composite throughput, averaged over each "
                 "grating's blazes", color=INK, fontsize=14, loc="left", pad=80)
    ax.annotate(
        f"composite = the MEAN of a grating's blaze curves at each wavelength, "
        f"then fitted with a cubic spline, knots every {KNOT_A:.0f} A\n"
        "each blaze fades in and out over the outer 15% of its own range, so "
        "one entering the average does not step the mean\n"
        "the whole fitted range is drawn; the steepest fit-end tails run off "
        "the bottom of the axis\n"
        "this is what the observed settings delivered, NOT what the grating can "
        "reach -- off-blaze configs pull it down; see kcrm_superblaze.png",
        xy=(0, 1.0), xycoords="axes fraction", textcoords="offset points",
        xytext=(0, 14), color=MUTED, fontsize=9.5, va="bottom", ha="left",
        annotation_clip=False)

    handles = [Line2D([], [], color=colour[g], lw=2.8, label=g) for g in ORDER]
    leg = ax.legend(handles=handles, frameon=False, fontsize=9.5,
                    loc="center left", bbox_to_anchor=(1.005, 0.5),
                    labelcolor=INK2)
    leg.set_title("grating", prop={"size": 9.5})
    leg.get_title().set_color(MUTED)
    if drawn_y:
        yy = np.concatenate(drawn_y)
        ax.set_ylim(ylow(yy), np.nanmax(yy) + 2.0)
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)

    with open(csv_path, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "wavelength_A", "composite_pct", "sd_pct",
                      "n_blazes", "n_effective"])
        for g, w, m, sd_, nn, ne in rows:
            wtr.writerow([g, f"{w:.0f}", f"{m:.3f}",
                          "" if not np.isfinite(sd_) else f"{sd_:.3f}", nn,
                          f"{ne:.2f}"])

    with open(os.path.join(ROOT, "composite_summary.csv"), "w",
              newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "nblazes", "composite_max_pct",
                      "at_wave_A", "n_min", "n_max", "small_drag_mean_pt",
                      "small_drag_worst_pt"])
        for r in summary:
            wtr.writerow([r["grating"], r["nblaze"], f"{r['peak']:.1f}",
                          f"{r['peak_wave']:.0f}", r["n_min"], r["n_max"],
                          "" if r["drag_mean"] is None else f"{r['drag_mean']:.2f}",
                          "" if r["drag_worst"] is None else f"{r['drag_worst']:.2f}"])
    return rows


def envelope(gg):
    """One grating's upper envelope: (grid, env, env_lim), or (None, None, None).

    The max over that grating's blaze means at each wavelength, then the spline.
    `env_lim` marks where the winning blaze is one whose every cube is flagged a
    lower limit -- a blaze like that is not a measurement of the grating, so the
    stretch it wins is drawn faint rather than read as the grating's efficiency.

    Factored out of superblaze() so anything comparing envelopes across
    reductions -- plot_notrim_compare.py -- is built from the same arithmetic
    and cannot drift from the figure it is compared against.
    """
    curves = blaze_means(gg)
    if not curves:
        return None, None, None
    lo = min(c[1].min() for c in curves)
    hi = max(c[1].max() for c in curves)
    grid = np.arange(lo, hi + 1.0, 1.0)
    stack, limit = [], []
    for deck, w, v, lim in curves:
        stack.append(np.interp(grid, w, v, left=np.nan, right=np.nan))
        limit.append(np.full(grid.shape, lim))
    stack, limit = np.vstack(stack), np.vstack(limit)
    ok = np.isfinite(stack).any(axis=0)
    env = np.full(grid.shape, np.nan)
    env_lim = np.zeros(grid.shape, dtype=bool)
    idx = np.nanargmax(np.where(np.isfinite(stack), stack, -np.inf), axis=0)
    env[ok] = stack[idx[ok], np.arange(grid.size)[ok]]
    env_lim[ok] = limit[idx[ok], np.arange(grid.size)[ok]]
    env[ok] = smooth(grid[ok], env[ok])
    return grid, env, env_lim


def superblaze(groups, path, title=None, extra_sub=None):
    """Per grating, the upper envelope of its blazes -- what it can reach.

    The envelope is the max over that grating's blaze means at each wavelength,
    i.e. what you would measure having tuned the grating to put its peak there.
    Where a Small-slicer blaze sets the maximum the line is dashed, because a
    narrower slicer reads ~6-8% low and that stretch understates the grating.
    """
    fig, ax = plt.subplots(figsize=(13.2, 6.6), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    colour = {g: HUES[i] for i, g in enumerate(ORDER)}
    rows, drawn_y = [], []

    for g in ORDER:
        gg = {k: v for k, v in groups.items() if k[0] == g}
        grid, env, env_lim = envelope(gg)
        if grid is None:
            continue
        ok = np.isfinite(env)

        # The slicer is NOT distinguished on this figure.  It is still part of
        # the blaze key, so no curve ever mixes deckers, and the per-blaze
        # figures and the tables carry it -- but splitting one grating's
        # envelope into two line styles asked the reader to hold a 7% caveat
        # while comparing seven gratings, and the caveat belongs in the text.
        # One stroke per grating, one alpha.  The `lower limit` flag is still
        # carried in blaze_summary.csv and blaze_peaks.csv, and it is still what
        # env_lim marks -- it is simply not drawn.  It applies to 69 A of RM1,
        # 2.6% of that curve, and a fade plus a legend row cost more of a
        # reader's attention than 69 A is worth on a figure whose job is the
        # SHAPE of seven curves.  NaNs outside the fitted range break the line
        # on their own.
        #
        # nm on x, fraction on y: the axes of the predicted-throughput figure
        # this is compared against.  Everything upstream -- peaks, tables, every
        # CSV -- stays in Angstrom and percent.
        ax.plot(grid / 10.0, env / 100.0, color=colour[g], lw=2.8, zorder=4,
                solid_capstyle="round")
        drawn_y.append(env[ok])

        for (gr, cw, deck), grp in sorted(gg.items(), key=lambda kv: kv[0][1]):
            kept = [c for c in grp if not c["drop"]]
            if not kept:
                continue
            gr2, mean2, _ = combine(grp)
            if gr2 is not None and len(kept) > 1:
                py, px, pl, ph = peak_stats(gr2, mean2)
                n = len(kept)
            else:
                c = max(kept, key=lambda c: c["peak"] if np.isfinite(c["peak"])
                        else -1)
                if not c["interior"].any():
                    continue
                py, px, pl, ph = peak_stats(c["wave"][c["interior"]],
                                            c["thru"][c["interior"]])
                n = 1
            if not np.isfinite(py):
                continue
            # No per-blaze marker: with the slicer distinction gone the figure
            # answers one question -- what does each grating reach -- and 30
            # dots answered a different one.  The peaks are still written to
            # blaze_peaks.csv and tabulated in GRATING_SUMMARY.md.
            rows.append(dict(grating=g, cenwave=cw, slicer=deck, ncubes=n,
                             peak=py, peak_wave=px, plat_lo=pl, plat_hi=ph,
                             limit=all(c["flag"] for c in kept)))

        anchor, dx, dy, ha = LABEL.get(g, (None, 0, 16, "center"))
        if anchor is None or not (grid[0] <= anchor <= grid[-1]):
            j = int(np.nanargmax(env))
        else:
            j = int(np.argmin(np.abs(grid - anchor)))
        ax.annotate(g, (grid[j] / 10.0, env[j] / 100.0),
                    textcoords="offset points",
                    xytext=(dx, dy), ha=ha, color=colour[g], fontsize=14,
                    fontweight="bold", zorder=8,
                    path_effects=[pe.withStroke(linewidth=3.5,
                                                foreground=SURFACE)])

    style_axes(ax, xlabel="Wavelength (nm)", ylabel="Throughput",
                label_fs=15, tick_fs=13)
    # The subtitle is three lines hung from the top of the axes; the title must
    # be padded clear of all of them or it lands on top of the first.
    ax.set_title(title or "KCRM red channel: what each grating reaches, "
                 "blaze by blaze",
                 color=INK, fontsize=14, loc="left", pad=80)
    sub = ("envelope = the best a grating measures at each wavelength, taken "
           "over its own central-wavelength settings, then fitted with a cubic "
        "spline, knots every 600 A\n"
           "an envelope falls at its own ends because only one setting reaches "
           "there, on its own falling edge -- not because the grating "
           "collapses\n"
        "nm and throughput-as-a-fraction, as the KCRM predicted-throughput "
        "figure plots them; y is cropped to 0.100-0.440, so the steepest "
        "fit-end tails run off the bottom")
    if extra_sub:
        sub += "\n" + extra_sub
    ax.annotate(sub, xy=(0, 1.0), xycoords="axes fraction",
                textcoords="offset points", xytext=(0, 14), color=MUTED,
                fontsize=9.5, va="bottom", ha="left", annotation_clip=False)

    handles = [Line2D([], [], color=colour[g], lw=2.8, label=g) for g in ORDER]
    leg = ax.legend(handles=handles, frameon=False, fontsize=9.5,
                    loc="center left", bbox_to_anchor=(1.005, 0.5),
                    labelcolor=INK2)
    leg.set_title("grating", prop={"size": 9.5})
    leg.get_title().set_color(MUTED)
    # Fixed to the reference figure's limits, not to what we drew: the point of
    # this variant is that the two can be laid side by side.  The terminal tails
    # that ylow() used to clip are visible again, which is the cost.
    ax.set_xlim(500, 1100)
    # y cropped to 0.100-0.440 rather than the reference figure's 0-0.5: the
    # curves live in 0.10-0.41, so the full range spent ~40% of the panel on
    # white space.  Terminal tails below 0.100 are clipped -- see the note in
    # the subtitle.
    ax.set_ylim(0.10, 0.44)
    ax.set_xticks([500, 700, 900, 1100])
    ax.set_yticks(np.arange(0.1, 0.41, 0.1))
    ax.set_yticklabels([f"{v:.3f}" for v in np.arange(0.1, 0.41, 0.1)])
    fig.savefig(path, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "blaze_plots"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cubes = collect()
    groups = collections.defaultdict(list)
    for c in cubes:
        groups[(c["grating"], c["cw"], c["decker"])].append(c)

    for g in ORDER:
        gg = {k: v for k, v in groups.items() if k[0] == g}
        if gg:
            grating_figure(g, gg,
                           os.path.join(args.outdir, f"{g.lower()}_blazes.png"))

    rows = superblaze(groups, os.path.join(ROOT, "kcrm_superblaze.png"))
    comp = composite(groups, os.path.join(ROOT, "kcrm_composite.png"),
                     os.path.join(ROOT, "composite_by_grating.csv"))
    trim_report(cubes, os.path.join(ROOT, "trim_summary.csv"))
    edge_residuals(os.path.join(ROOT, "edge_check.csv"))

    # Per-cube CSV: the table in GRATING_SUMMARY.md is generated from this, so
    # the document and the figures cannot drift apart.
    csv_path = os.path.join(ROOT, "blaze_summary.csv")
    with open(csv_path, "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "cenwave", "slicer", "binning", "night",
                      "setup", "star", "nframes", "exptime_s", "airmass",
                      "fit_lo_A", "fit_hi_A", "peak_pct", "peak_wave_A",
                      "interior_mean_pct", "status"])
        for key in sorted(groups, key=lambda k: (ORDER.index(k[0]), k[1])):
            for c in sorted(groups[key], key=lambda c: c["night"]):
                st = "excluded" if c["drop"] else ("lower limit" if c["flag"]
                                                   else "ok")
                wtr.writerow([c["grating"], c["cw"], c["decker"], c["binning"],
                              c["night"], c["setup"], c["star"], c["nframes"],
                              f"{c['exptime']:.0f}", f"{c['airmass']:.2f}",
                              f"{c['lo']:.0f}", f"{c['hi']:.0f}",
                              f"{c['peak']:.1f}", f"{c['peak_wave']:.0f}",
                              f"{c['mean']:.1f}", st])

    with open(os.path.join(ROOT, "blaze_peaks.csv"), "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["grating", "cenwave", "slicer", "ncubes",
                      "combined_peak_pct", "peak_wave_A", "plateau_lo_A",
                      "plateau_hi_A", "status"])
        for r in rows:
            wtr.writerow([r["grating"], r["cenwave"], r["slicer"], r["ncubes"],
                          f"{r['peak']:.1f}", f"{r['peak_wave']:.0f}",
                          f"{r['plat_lo']:.0f}", f"{r['plat_hi']:.0f}",
                          "lower limit" if r["limit"] else "ok"])

    print(f"{len(cubes)} cubes, {len(groups)} blazes, "
          f"{len({r[0] for r in comp})} composites")
    for r in rows:
        print(f"  {r['grating']:4s} {r['cenwave']:5d} {r['slicer']:6s} "
              f"n={r['ncubes']}  peak {r['peak']:5.1f}% @ {r['peak_wave']:.0f} A"
              f"  (plateau {r['plat_lo']:.0f}-{r['plat_hi']:.0f} A)")


if __name__ == "__main__":
    main()
