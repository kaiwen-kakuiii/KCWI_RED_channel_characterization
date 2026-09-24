#!/usr/bin/env python
"""Draw where every config's spectrum was cut, and write the ledger behind it.

Each reduced cube is narrowed twice before any throughput number is read off it:

  delivered   the extracted spectrum's own extent (spec1d BOX_WAVE), which is the
              datacube's wavelength axis -- extraction cuts spatially, never in
              wavelength
  after trim  SENS_WAVE's non-zero extent, i.e. what pypeit_sensfunc was actually
              handed after measure_trim.py's counts + full-slice-coverage cut
  fit         SENS_ZEROPOINT_FIT_GPM's extent, what survived the zeropoint fit

The cut is read off those three products, not off a parameter file, so what is
plotted is what was applied.  The declared value is carried alongside as a
cross-check, read from the per-product `.sens` file -- NOT from `sensfunc.par`,
which pypeit rewrites on every run and which therefore holds only the last
variant fitted in that directory.  Any config where declared and applied
disagree by more than 3 px is flagged.

Outputs:
    kcrm_trim_map.png   per config: delivered vs kept in Angstrom, and the cut in pixels
    trim_detail.csv     the same, per config, with both ends separated

Usage:
    python plot_trim_map.py [--outdir <dir>]
"""
import argparse
import csv
import glob
import os
import re

import numpy as np
from astropy.io import fits

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import Patch                  # noqa: E402

import plot_blaze_summary as P

ROOT = P.ROOT
ORDER, HUES = P.ORDER, P.HUES
INK, INK2, MUTED = P.INK, P.INK2, P.MUTED
GRID, SPINE, SURFACE = P.GRID, P.SPINE, P.SURFACE

PAR_TOL = 3          # px; above this the declared and applied cuts disagree


def declared_trim(sdir, tag, suffix):
    """`trim_std_pixs = blue, red` as declared for THIS product, or (-1, -1).

    The record is the `.sens` parameter file the product was fitted with, named
    per product and per variant.  `sensfunc.par` in the same directory is not
    it: pypeit rewrites that file on every run, so where a directory holds
    several products or several polynomial orders it keeps only the last one --
    RL 2024-12-28 B's sensfunc.par still carries the p15 run's `10, 146` while
    its p15cov product was fitted with `261, 259`.  Matched by name first, then
    a lone file, never an ambiguous one.
    """
    cand = os.path.join(sdir, f"{tag}_{suffix}.sens")
    if not os.path.exists(cand):
        hits = [f for f in glob.glob(os.path.join(sdir, f"*_{suffix}.sens"))]
        if len(hits) != 1:
            return -1, -1
        cand = hits[0]
    m = re.search(r"trim_std_pixs\s*=\s*([0-9]+)\s*,\s*([0-9]+)",
                  open(cand).read())
    return (int(m.group(1)), int(m.group(2))) if m else (-1, -1)


def ledger():
    """One row per surviving config: what it delivered and what was cut off it."""
    rows = []
    for g in ORDER:
        pat = os.path.join(ROOT, f"{g} pypeit run", "*", "pypeit_run",
                           "keck_kcrm_*", f"*_sens_{P.SUFFIX[g]}.fits")
        for sens in sorted(glob.glob(pat)):
            if "order_study" in sens:
                continue
            sdir = os.path.dirname(sens)
            tag = os.path.basename(sens).split("_sens_")[0]
            c = P.read_cube(sens, g)
            if c is None or c["drop"] or not np.isfinite(c["deliv_lo"]):
                continue
            s1 = os.path.join(sdir, f"spec1d_{tag}.fits")
            w1 = np.asarray(fits.getdata(s1, 1)["BOX_WAVE"]).flatten()
            w1 = np.sort(w1[w1 > 1])
            n = w1.size
            # The applied cut, expressed in pixels of the delivered spectrum.
            b_px = int(np.searchsorted(w1, c["sens_lo"]))
            r_px = int(n - np.searchsorted(w1, c["sens_hi"], side="right"))
            pb, pr = declared_trim(sdir, tag, P.SUFFIX[g])
            rows.append(dict(
                grating=g, cw=c["cw"], slicer=c["decker"], night=c["night"],
                setup=c["setup"], star=c["star"], npix=n,
                disp=float(np.median(np.diff(w1))),
                blue_px=b_px, red_px=r_px,
                blue_A=c["sens_lo"] - c["deliv_lo"],
                red_A=c["deliv_hi"] - c["sens_hi"],
                fitcut_blue_A=c["lo"] - c["sens_lo"],
                fitcut_red_A=c["sens_hi"] - c["hi"],
                par_blue=pb, par_red=pr,
                par_ok=abs(pb - b_px) <= PAR_TOL and abs(pr - r_px) <= PAR_TOL,
                deliv_lo=c["deliv_lo"], deliv_hi=c["deliv_hi"],
                sens_lo=c["sens_lo"], sens_hi=c["sens_hi"],
                fit_lo=c["lo"], fit_hi=c["hi"], tag=tag))
    rows.sort(key=lambda r: (ORDER.index(r["grating"]), r["cw"], r["night"]))
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def figure(rows, path):
    """Left: the cut in wavelength.  Right: the same cut in pixels."""
    hue = {g: HUES[i] for i, g in enumerate(ORDER)}

    # One y slot per config, with a blank slot between gratings.
    y, lab, seen = [], [], None
    slot = 0
    for r in rows:
        if seen is not None and r["grating"] != seen:
            slot += 1.2
        seen = r["grating"]
        y.append(slot)
        lab.append(f"{r['cw']} {r['slicer'][0]}  {r['night']}")
        slot += 1
    y = np.array(y)

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(15.5, 13.0), sharey=True,
        gridspec_kw=dict(width_ratios=[3.05, 1], wspace=0.045))
    fig.patch.set_facecolor(SURFACE)

    for a in (ax, bx):
        a.set_facecolor(SURFACE)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_color(SPINE)
        a.tick_params(colors=INK2, labelsize=8, length=3, color=SPINE)
        a.set_axisbelow(True)

    h = 0.62
    for r, yy in zip(rows, y):
        c = hue[r["grating"]]
        # Delivered: the whole spectrum, recessive.
        ax.barh(yy, r["deliv_hi"] - r["deliv_lo"], left=r["deliv_lo"], height=h,
                color=GRID, edgecolor=SPINE, linewidth=0.4, zorder=2)
        # Kept: what the throughput curve is actually read from.
        ax.barh(yy, r["fit_hi"] - r["fit_lo"], left=r["fit_lo"], height=h,
                color=c, edgecolor=SURFACE, linewidth=0.8, zorder=3)
        # How much each end gave up, in Angstrom.
        ax.text(r["deliv_lo"] - 40, yy, f"{r['blue_A']:.0f}", ha="right",
                va="center", fontsize=6.4, color=MUTED, zorder=4)
        ax.text(r["deliv_hi"] + 40, yy, f"{r['red_A']:.0f}", ha="left",
                va="center", fontsize=6.4, color=MUTED, zorder=4)

        # The same cut as a fraction of the spectrum, blue end then red end.
        # Raw pixel counts are not comparable across configs: the Small-slicer
        # nights are unbinned, ~4600 px against ~2200 for the binned ones.
        fb = 100.0 * r["blue_px"] / r["npix"]
        fr = 100.0 * r["red_px"] / r["npix"]
        bx.barh(yy, fb, height=h, color=c, zorder=3)
        bx.barh(yy, fr, left=fb, height=h, color=c, alpha=0.42,
                edgecolor=SURFACE, linewidth=0.8, zorder=3)
        bx.text(fb + fr + 0.5, yy, f"{r['blue_px']}+{r['red_px']} px",
                ha="left", va="center", fontsize=6.4, color=MUTED, zorder=4)

    # Grating name beside its block, in its own hue.
    for g in ORDER:
        idx = [i for i, r in enumerate(rows) if r["grating"] == g]
        yc = 0.5 * (y[idx[0]] + y[idx[-1]])
        cut = np.mean([rows[i]["blue_A"] + rows[i]["red_A"] for i in idx])
        px = np.mean([100.0 * (rows[i]["blue_px"] + rows[i]["red_px"])
                      / rows[i]["npix"] for i in idx])
        ax.text(4460, yc, g, ha="left", va="center", fontsize=13,
                fontweight="bold", color=hue[g])
        ax.text(4460, yc + 0.78, f"−{cut:.0f} A", ha="left", va="center",
                fontsize=7.5, color=MUTED)
        ax.text(4460, yc + 1.62, f"{px:.0f}% cut", ha="left", va="center",
                fontsize=7.5, color=MUTED)

    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=6.6, color=INK2, fontfamily="monospace")
    ax.invert_yaxis()
    ax.set_xlim(4400, 10950)
    ax.set_xlabel("wavelength (A)", fontsize=9.5, color=INK2)
    ax.xaxis.grid(True, color=GRID, linewidth=0.7)
    ax.set_title("What each config delivered, and what survived to the fit",
                 fontsize=12, color=INK, loc="left", pad=26)
    ax.text(0, 1.012, "grey = delivered by the reduction   ·   colour = kept "
                      "(the range every throughput number is read from)   ·   "
                      "numbers outside each bar = A cut off that end",
            transform=ax.transAxes, fontsize=8.2, color=MUTED)

    bx.set_xlim(0, 45)
    bx.set_xlabel("% of the spectrum cut  (solid = blue end, pale = red end)",
                  fontsize=9.5, color=INK2)
    bx.xaxis.grid(True, color=GRID, linewidth=0.7)
    bx.set_title("The same cut, in % of the spectrum",
                 fontsize=12, color=INK, loc="left", pad=26)
    bx.text(0, 1.012, "11–25% everywhere: one operation, seven gratings",
            transform=bx.transAxes, fontsize=8.2, color=MUTED)

    fig.suptitle("KCRM red channel — where each config's spectrum was trimmed",
                 fontsize=15.5, color=INK, x=0.077, ha="left", y=0.986)
    fig.text(0.077, 0.9585,
             "Cut read off the products: spec1d BOX_WAVE (delivered) vs "
             "SENS_ZEROPOINT_FIT_GPM (kept).  measure_trim.py does almost all of it\n"
             "— a counts threshold AND full 24-slice coverage — and the zeropoint "
             "fit removes 0–13 A more.",
             fontsize=8.6, color=MUTED, ha="left", va="top")

    fig.subplots_adjust(left=0.152, right=0.972, top=0.902, bottom=0.048)
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=ROOT)
    a = ap.parse_args()
    rows = ledger()
    write_csv(rows, os.path.join(a.outdir, "trim_detail.csv"))
    figure(rows, os.path.join(a.outdir, "kcrm_trim_map.png"))
    bad = [f"{r['grating']} {r['tag']}" for r in rows if not r["par_ok"]]
    print(f"{len(rows)} configs -> trim_detail.csv, kcrm_trim_map.png")
    if bad:
        print("declared .sens trim disagrees with the product: " + ", ".join(bad))


if __name__ == "__main__":
    main()
