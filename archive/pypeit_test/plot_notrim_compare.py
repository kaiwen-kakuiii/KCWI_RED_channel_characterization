#!/usr/bin/env python
"""The superblaze WITHOUT measure_trim.py's cut, against the one with it.

The comparison exists because the two figures are not two views of one dataset:
the *_notrim products are a separate sensfunc fit (run_notrim_sens.py) over the
reduction's full delivered range.  The only cut left in them is pypeit's own --
the sigma-rejection inside the zeropoint fit, recorded in
SENS_ZEROPOINT_FIT_GPM -- so the difference between the two curves is exactly
what our trim removes, plus the change in the polynomial that follows from
fitting over a wider range.

That second term is not a detail.  The polynomial is fitted over the whole
input, so restoring a biased zone does not merely append it to the end of the
curve: it moves the fit everywhere, including the interior that both variants
share.  Any interior difference here is that effect, not a measurement.

Three outputs:

  kcrm_superblaze_notrim.png    the no-trim envelope alone, drawn like the real
                                one so the two files can be flipped between
  kcrm_superblaze_compare.png   per grating, both envelopes on one axis
  notrim_compare.csv            span, peak and peak wavelength, both variants

Usage:
    python plot_notrim_compare.py
"""
import collections
import csv
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402

import plot_blaze_summary as P

ROOT, ORDER, HUES = P.ROOT, P.ORDER, P.HUES
INK, INK2, MUTED = P.INK, P.INK2, P.MUTED
GRID, SPINE, SURFACE = P.GRID, P.SPINE, P.SURFACE

NOTRIM = {g: s.replace("cov", "notrim") for g, s in P.SUFFIX.items()}


def groups_for(suffix):
    """Blaze groups built from one set of sensfunc products.

    P.collect() reads P.SUFFIX at call time, so swapping it is enough -- and it
    guarantees both variants go through identical code, which is the only way
    the difference between the two figures means what it claims to.
    """
    keep = dict(P.SUFFIX)
    try:
        P.SUFFIX.update(suffix)
        cubes = P.collect()
    finally:
        P.SUFFIX.clear()
        P.SUFFIX.update(keep)
    groups = collections.defaultdict(list)
    for c in cubes:
        groups[(c["grating"], c["cw"], c["decker"])].append(c)
    return cubes, groups


def envelopes(groups):
    """{grating: (grid, env)} over the gratings that have any blaze."""
    out = {}
    for g in ORDER:
        gg = {k: v for k, v in groups.items() if k[0] == g}
        grid, env, _ = P.envelope(gg)
        if grid is not None:
            out[g] = (grid, env)
    return out


def compare_figure(cov, notrim, path):
    """One small multiple per grating: both envelopes, same axes, same hue."""
    hue = {g: HUES[i] for i, g in enumerate(ORDER)}
    fig, axes = plt.subplots(4, 2, figsize=(13.4, 14.6), dpi=170)
    fig.patch.set_facecolor(SURFACE)

    for ax, g in zip(axes.ravel(), ORDER):
        ax.set_facecolor(SURFACE)
        c = hue[g]
        gc, ec = cov.get(g, (None, None))
        gn, en = notrim.get(g, (None, None))
        if gn is not None:
            ax.plot(gn, en, color=c, lw=2.2, ls=(0, (5, 2.4)), zorder=3,
                    label="no trim (pypeit's fit rejection only)")
        if gc is not None:
            ax.plot(gc, ec, color=c, lw=2.8, zorder=4, solid_capstyle="round",
                    label="as published (measure_trim applied)")
        # Shade what the trim removed, at each end.
        if gc is not None and gn is not None:
            for a, b in ((gn[0], gc[0]), (gc[-1], gn[-1])):
                if b > a:
                    ax.axvspan(a, b, color=c, alpha=0.10, lw=0, zorder=1)
        ax.set_title(g, color=c, fontsize=12.5, fontweight="bold", loc="left",
                     pad=6)
        if gc is not None and gn is not None:
            # Usually the no-trim fit is wider.  It is not always: handed the
            # full range, pypeit's own sigma-rejection can throw away MORE than
            # our trim did, which is worth saying in words rather than as a
            # negative "restored".
            dd = (gn[-1] - gn[0]) - (gc[-1] - gc[0])
            d = (f"+{dd:.0f} A restored" if dd >= 0 else
                 f"{-dd:.0f} A NARROWER — pypeit rejected more than we trimmed")
            pk = f"peak {np.nanmax(ec):.1f}% → {np.nanmax(en):.1f}%"
            ax.set_title(f"{d}   ·   {pk}", color=MUTED, fontsize=9,
                         loc="right", pad=8)
        ax.grid(color=GRID, linewidth=0.7)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(SPINE)
        ax.tick_params(colors=INK2, labelsize=8.5, length=3, color=SPINE)
        yy = np.concatenate([v for v in (ec, en) if v is not None])
        ax.set_ylim(P.ylow(yy), np.nanmax(yy) + 2.0)
        ax.set_xlabel("wavelength (A)", fontsize=9, color=INK2)
        ax.set_ylabel("throughput (%)", fontsize=9, color=INK2)

    last = axes.ravel()[len(ORDER)]
    last.axis("off")
    last.legend(handles=[
        Line2D([], [], color=INK2, lw=2.8, label="as published — "
               "measure_trim.py applied"),
        Line2D([], [], color=INK2, lw=2.2, ls=(0, (5, 2.4)),
               label="no trim — pypeit's fit rejection only"),
        Line2D([], [], color=INK2, lw=8, alpha=0.10,
               label="the range the trim removes")],
        frameon=False, fontsize=10, loc="upper left", labelcolor=INK2,
        bbox_to_anchor=(0.02, 0.86))
    last.text(0.02, 0.94, "what the two lines are", transform=last.transAxes,
              fontsize=10.5, color=INK, fontweight="bold")
    last.text(0.02, 0.36,
              "Both are the same envelope: per grating, the max over its\n"
              "central-wavelength settings, spline-smoothed, knots every 600 A.\n"
              "The dashed line is a separate sensfunc fit over the full\n"
              "delivered range, so it also differs in the interior — a wider\n"
              "fit moves the polynomial everywhere, not only at the ends.",
              transform=last.transAxes, fontsize=8.8, color=MUTED, va="top")

    fig.suptitle("KCRM red channel — the superblaze with and without our trim",
                 fontsize=15, color=INK, x=0.055, ha="left", y=0.988)
    fig.text(0.055, 0.968,
             "The restored range is the partial-slice-coverage zone: near a "
             "cube's ends only some of the 24 slices contribute, and the coadd "
             "fills the rest by smearing\nneighbouring slices' flux. Counts "
             "there look normal and the zeropoint is biased up to ~0.2 mag. "
             "That is why the solid line is the published one.",
             fontsize=9, color=MUTED, ha="left", va="top")
    fig.subplots_adjust(left=0.058, right=0.982, top=0.925, bottom=0.042,
                        hspace=0.36, wspace=0.19)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def table(cov, notrim, path):
    rows = []
    for g in ORDER:
        if g not in cov or g not in notrim:
            continue
        gc, ec = cov[g]
        gn, en = notrim[g]
        rows.append(dict(
            grating=g,
            cov_lo=gc[0], cov_hi=gc[-1], cov_span=gc[-1] - gc[0],
            notrim_lo=gn[0], notrim_hi=gn[-1], notrim_span=gn[-1] - gn[0],
            restored=(gn[-1] - gn[0]) - (gc[-1] - gc[0]),
            cov_peak=np.nanmax(ec), cov_peak_wave=gc[np.nanargmax(ec)],
            notrim_peak=np.nanmax(en), notrim_peak_wave=gn[np.nanargmax(en)],
            # Interior = the range both cover.  A difference here is the fit
            # moving under a wider input, not a different measurement.
            interior_mean_delta=float(np.nanmean(
                np.interp(gc, gn, en) - ec)),
            interior_max_delta=float(np.nanmax(
                np.interp(gc, gn, en) - ec))))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.1f}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    return rows


def main():
    _, gc = groups_for(P.SUFFIX)
    _, gn = groups_for(NOTRIM)
    cov, notrim = envelopes(gc), envelopes(gn)

    P.superblaze(
        gn, os.path.join(ROOT, "kcrm_superblaze_notrim.png"),
        title="KCRM red channel: what each grating reaches — NO TRIM "
              "(comparison only)",
        extra_sub="measure_trim.py's cut removed: the full delivered range is "
                  "fitted, so the ends sit in the partial-slice-coverage zone "
                  "and are biased up to ~0.2 mag")
    compare_figure(cov, notrim,
                   os.path.join(ROOT, "kcrm_superblaze_compare.png"))
    rows = table(cov, notrim, os.path.join(ROOT, "notrim_compare.csv"))

    print(f"{'g':4} {'span cov':>9} {'span raw':>9} {'restored':>9} "
          f"{'peak cov':>9} {'peak raw':>9} {'interior d':>11} {'max d':>7}")
    for r in rows:
        print(f"{r['grating']:4} {r['cov_span']:9.0f} {r['notrim_span']:9.0f} "
              f"{r['restored']:9.0f} {r['cov_peak']:8.1f}% {r['notrim_peak']:8.1f}% "
              f"{r['interior_mean_delta']:+11.2f} {r['interior_max_delta']:+7.2f}")


if __name__ == "__main__":
    main()
