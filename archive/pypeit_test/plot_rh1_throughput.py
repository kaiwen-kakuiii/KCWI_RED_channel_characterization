#!/usr/bin/env python
"""Plot every throughput curve of one grating on one axis, RL underneath for scale.

Each cube is drawn only across its own fitted range (SENS_ZEROPOINT_FIT_GPM);
outside it the zeropoint polynomial extrapolates and has been seen to exceed
100%.  The outer 150 A of each range is drawn faint, because the fit is worth
~3% in the interior and up to ~12% within 150 A of an end.

Colour encodes CENTRAL WAVELENGTH, not night.  Changing central wavelength
rotates the grating and moves the blaze peak, so two setups sample a different
efficiency curve even at the same wavelength -- cenwave is the variable that
makes two curves comparable or not.  Airmass, the other thing that separates
cubes of one setup, is carried by line style rather than a second hue, so the
two never compete for the same channel.

Defaults to RH1; pass --grating RH2 for the RH2 run.  Slicer is deliberately NOT
encoded: in the RH2 set every Medium config sits at cenwave 6750 and every Large
at 6900/7100, so slicer and central wavelength are perfectly confounded and a
slicer colouring would assert a distinction the data cannot support.

Nights in rh1_throughput_table.DROP are not drawn -- one source of truth, so the
figure and the table can never disagree about what was excluded.  --keep-dropped
draws them anyway, in grey, to see what was removed.

Usage:
    python plot_rh1_throughput.py [--grating RH1|RH2] [-o out.png] [--keep-dropped]
"""
import argparse
import collections
import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RL_CSV = os.path.join(ROOT, "throughput_composite_p15cov.csv")

sys.path.insert(0, HERE)
from rh1_throughput_table import DROP                  # noqa: E402

# Validated categorical palette (dataviz reference instance, all checks pass).
# Assigned in fixed order to central wavelength, ascending -- never cycled.
#
# A SEVENTH slot (violet) was appended for RM1, which has seven central
# wavelengths where every earlier grating had at most four.  Appending rather
# than re-ordering is deliberate: hue must follow the entity, so the first six
# assignments -- and therefore every RH1, RH2 and RH3 figure already on disk --
# are byte-identical before and after.  Re-validated as a 7-slot set on the
# light surface #fcfcfb: worst adjacent CVD dE 9.1 (protan), worst adjacent
# normal-vision dE 19.6, both clear.  The contrast check WARNs on aqua, yellow
# and magenta, which obligates relief -- the legend below and the results table
# in each grating's .md are it.
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
        "#4a3aa7"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
SURFACE = "#fcfcfb"
EDGE = 150.0


def collect(runs, pattern="*_sens_*.fits"):
    out = []
    for sens in sorted(glob.glob(os.path.join(runs, "*", "pypeit_run",
                                              "keck_kcrm_*", pattern))):
        tag = os.path.basename(sens).split("_sens_")[0]
        sdir = os.path.dirname(sens)
        with fits.open(sens) as h:
            t = Table(h["SENS"].data)
            wave = np.asarray(h["WAVE"].data).flatten()
            thru = np.asarray(h["THROUGHPUT"].data).flatten() * 100.0
        gpm = np.asarray(t["SENS_ZEROPOINT_FIT_GPM"][0], dtype=bool)
        sw = np.asarray(t["SENS_WAVE"][0]).flatten()
        lo, hi = sw[gpm].min(), sw[gpm].max()
        keep = (wave >= lo) & (wave <= hi) & np.isfinite(thru) & (thru > 0)
        if not keep.any():
            continue

        cw, am = np.nan, np.nan
        pf = glob.glob(os.path.join(sdir, "*.pypeit"))
        if pf:
            m = re.search(r"cenwave:\s*([0-9.]+)", open(pf[0]).read())
            if m:
                # Round: the same nominal grating setting is written 6519.939 on
                # one night and 6519.996 on another, which would otherwise split
                # one central wavelength into two colours.
                cw = round(float(m.group(1)))
        cube = os.path.join(sdir, f"{tag}.fits")
        if os.path.exists(cube):
            am = fits.getheader(cube).get("AIRMASS", np.nan)
        out.append(dict(tag=tag, wave=wave[keep], thru=thru[keep],
                        lo=lo, hi=hi, cw=cw, am=am,
                        drop=next((n for n in DROP if n in tag), None)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sens-glob", default="*_sens_*.fits",
                    help="which sensfunc files to plot.  Needed when several "
                         "polynomial orders coexist on disk: RH3 keeps p3/p5/p7/"
                         "p9/p15 from its order study, and the default glob would "
                         "draw all of them as separate curves")
    ap.add_argument("-g", "--grating", default="RH1",
                    help="which run directory to read, e.g. RH1 or RH2")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--keep-dropped", action="store_true",
                    help="also draw the excluded nights, in grey")
    args = ap.parse_args()
    runs = os.path.join(ROOT, f"{args.grating} pypeit run")
    if args.out is None:
        args.out = os.path.join(ROOT,
                                f"{args.grating.lower()}_throughput_curves.png")

    cubes = collect(runs, args.sens_glob)
    if not cubes:
        sys.exit("no sensfunc files found")

    # Assign hues over EVERY cenwave present, including the dropped nights, then
    # filter.  Colour must follow the entity, not its rank: 2024-03-13 is the
    # only cube at 6150 A, so keying off the surviving list would shift every
    # cenwave above it one hue along and repaint curves that did not change.
    cws = sorted({c["cw"] for c in cubes if np.isfinite(c["cw"])})
    colour = {cw: HUES[i % len(HUES)] for i, cw in enumerate(cws)}

    dropped = [c for c in cubes if c["drop"]]
    if not args.keep_dropped:
        cubes = [c for c in cubes if not c["drop"]]

    fig, ax = plt.subplots(figsize=(12.4, 6.2), dpi=170)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # RL for scale, recessive: a different grating over a different range, shown
    # so the RH1 numbers are read against something, not mistaken for it.
    if os.path.exists(RL_CSV):
        import csv
        w, v, s = [], [], []
        with open(RL_CSV) as fh:
            for r in csv.DictReader(fh):
                w.append(float(r["wavelength_A"]))
                v.append(float(r["throughput_smooth"]) * 100)
                s.append(float(r["throughput_std_smooth"]) * 100)
        w, v, s = np.array(w), np.array(v), np.array(s)
        ax.fill_between(w, v - s, v + s, color=MUTED, alpha=0.18, lw=0, zorder=1)
        ax.plot(w, v, color=MUTED, lw=2, zorder=2)
        rl_w, rl_v = w, v

    byc = collections.defaultdict(list)
    for c in cubes:
        byc[c["cw"]].append(c)

    for cw, group in sorted(byc.items()):
        group.sort(key=lambda c: c["am"] if np.isfinite(c["am"]) else 0)
        for k, c in enumerate(group):
            col = MUTED if c["drop"] else colour.get(cw, MUTED)
            style = "-" if k == 0 else "--"
            m = (c["wave"] > c["lo"] + EDGE) & (c["wave"] < c["hi"] - EDGE)
            ax.plot(c["wave"], c["thru"], style, color=col, lw=2,
                    alpha=0.30, zorder=3)                    # edges, faint
            if m.any():
                ax.plot(c["wave"][m], c["thru"][m], style, color=col, lw=2,
                        alpha=1.0, zorder=4)                 # interior, solid

    # No direct labels: there are six or seven central wavelengths and the
    # method allows direct labelling only up to four -- past that they collide
    # with the curves they name.  The legend carries identity instead, and the
    # palette passed the normal-vision and CVD separation checks, so hue alone
    # is legible here.

    ax.set_xlabel("wavelength (A)", color=INK2, fontsize=10)
    ax.set_ylabel("throughput (%)", color=INK2, fontsize=10)
    # Title pad has to clear the subtitle, which is two lines whenever nights
    # are excluded.  Anchor the subtitle from its TOP and hang it downward, so
    # adding a line grows it toward the axes instead of up through the title.
    nsub = 2 if (dropped and not args.keep_dropped) else 1
    ax.set_title(f"KCRM {args.grating} throughput, one curve per visit",
                 color=INK, fontsize=13, loc="left", pad=16 + 15 * nsub)
    sub = ("colour = grating central wavelength   ·   solid = lower airmass, "
           "dashed = higher   ·   faint = within 150 A of a fit edge")
    if dropped and not args.keep_dropped:
        sub += ("\n" + ", ".join(sorted(DROP)) +
                " excluded: non-photometric, achromatic light loss at identical "
                "instrument config")
    ax.annotate(sub, xy=(0, 1.0), xycoords="axes fraction",
                xytext=(0, 14 + 15 * nsub), textcoords="offset points",
                color=MUTED, fontsize=9, va="top", ha="left",
                annotation_clip=False)

    lo = min(c["lo"] for c in cubes) - 60
    hi = max(c["hi"] for c in cubes) + 220
    ax.set_xlim(lo, hi)
    # y-range from the INTERIOR of each fit only: the faint edge tails plunge
    # toward zero and would otherwise squash every curve into the top third.
    iv = np.concatenate([c["thru"][(c["wave"] > c["lo"] + EDGE) &
                                   (c["wave"] < c["hi"] - EDGE)]
                         for c in cubes if ((c["wave"] > c["lo"] + EDGE) &
                                            (c["wave"] < c["hi"] - EDGE)).any()])
    ax.set_ylim(max(0, iv.min() - 6), iv.max() + 4)

    # Anchor the RL label inside the x-range actually being drawn.
    if "rl_w" in dir():
        xa = lo + 0.06 * (hi - lo)
        if rl_w.min() <= xa <= rl_w.max():
            ja = int(np.argmin(np.abs(rl_w - xa)))
            ax.annotate("RL grating (previous work)", (rl_w[ja], rl_v[ja]),
                        textcoords="offset points", xytext=(4, -22),
                        color=INK2, fontsize=9, ha="left")

    ax.grid(True, color="#e6e5e0", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#d8d7d2")
    ax.tick_params(colors=INK2, labelsize=9)

    # Legend from the cenwaves actually drawn, not from `cws`: dropping
    # 2024-03-13 removes 6150 A entirely, and a legend entry with no curve on
    # the axis sends the reader hunting for something that is not there.
    drawn = sorted({c["cw"] for c in cubes if np.isfinite(c["cw"])})
    handles = [Line2D([], [], color=colour[cw], lw=2, label=f"{cw:.0f} A")
               for cw in drawn]
    handles += [Line2D([], [], color=INK2, lw=2, ls="-", label="lower airmass"),
                Line2D([], [], color=INK2, lw=2, ls="--", label="higher airmass")]
    if dropped and args.keep_dropped:
        handles += [Line2D([], [], color=MUTED, lw=2, label="excluded night")]
    leg = ax.legend(handles=handles, frameon=False, fontsize=9, ncol=1,
                    loc="center left", bbox_to_anchor=(1.005, 0.5),
                    labelcolor=INK2)
    leg.set_title("central wavelength", prop={"size": 9})
    leg.get_title().set_color(MUTED)

    fig.tight_layout()
    # bbox_inches="tight": the subtitle is an annotation outside the axes and
    # tight_layout does not reserve room for it, so it would be cropped.
    fig.savefig(args.out, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {args.out}  ({len(cubes)} cubes, "
          f"{len(drawn)} central wavelengths"
          + (f", {len(dropped)} excluded" if dropped and not args.keep_dropped
             else "") + ")")

    print(f"\n{'cube':30s} {'cenw':>6s} {'airm':>5s} {'range (A)':>16s} {'interior mean':>14s}")
    for c in sorted(cubes, key=lambda c: (c["cw"], c["am"])):
        m = (c["wave"] > c["lo"] + EDGE) & (c["wave"] < c["hi"] - EDGE)
        v = f"{c['thru'][m].mean():13.1f}%" if m.any() else "             -"
        print(f"{c['tag']:30s} {c['cw']:6.0f} {c['am']:5.2f} "
              f"{c['lo']:7.1f}-{c['hi']:7.1f} {v}")


if __name__ == "__main__":
    main()
