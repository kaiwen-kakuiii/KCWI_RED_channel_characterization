#!/usr/bin/env python
"""Overlay every config's sensitivity function on one figure, grouped by slicer.

The RL equivalent (all_configs_sens_overlay_p15cov.png) was made ad hoc with no
script; this is that plot, generalised, so RL and RH2 are drawn the same way.

Two panels sharing the wavelength axis:
    top     zeropoint: measured points inside the fit range (faint) + the fit
    bottom  throughput, drawn only inside the fitted range

**Colour encodes the slicer, not the config.**  On RH2 the slicer is the axis the
result splits along (Large ~30%, Medium ~24%), and six arbitrary hues would hide
exactly the thing worth seeing.  Configs within a slicer are separated by line
style instead.  Everything is masked to SENS_ZEROPOINT_FIT_GPM: outside it the
polynomial extrapolates and can be wrong by magnitudes.

Usage:
    python plot_sens_overlay.py "<runroot>" [suffix] [-o out.png]
        runroot   e.g. "RH2 pypeit run"
        suffix    sens tag, default p15cov  (matches *_sens_IR_p15cov.fits)
"""
import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Categorical slots 1 and 2 of the validated default palette.  Checked with
# scripts/validate_palette.js --mode light --pairs all: CVD dE 24.7, normal 33.6,
# both >= 3:1 on the surface.  Do not swap for arbitrary hues.
SLICER_COLOR = {"Large": "#2a78d6", "Medium": "#eb6834", "Small": "#1baf7a"}
STYLES = ["-", "--", ":", "-."]
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"


def load(path):
    with fits.open(path) as h:
        d = h["SENS"].data
        w = np.asarray(d["SENS_WAVE"]).ravel()
        zp = np.asarray(d["SENS_ZEROPOINT"]).ravel()
        zf = np.asarray(d["SENS_ZEROPOINT_FIT"]).ravel()
        gpm = np.asarray(d["SENS_ZEROPOINT_FIT_GPM"]).ravel().astype(bool)
        wv = np.asarray(h["WAVE"].data).ravel()
        th = np.asarray(h["THROUGHPUT"].data).ravel()
    ok = gpm & np.isfinite(w) & (w > 0) & np.isfinite(zf)
    if not ok.any():
        return None
    lo, hi = w[ok].min(), w[ok].max()
    tm = np.isfinite(wv) & np.isfinite(th) & (wv >= lo) & (wv <= hi) & (th > 0)
    return dict(w=w[ok], zp=zp[ok], zf=zf[ok], lo=lo, hi=hi,
                tw=wv[tm], th=100 * th[tm])


def main():
    a = [x for x in sys.argv[1:]]
    out = None
    if "-o" in a:
        i = a.index("-o"); out = a[i + 1]; del a[i:i + 2]
    if not a:
        sys.exit(__doc__)
    runroot = a[0]
    suffix = a[1] if len(a) > 1 else "p15cov"

    pat = os.path.join(runroot, "*", "pypeit_run", "keck_kcrm_*",
                       f"*_sens_IR_{suffix}.fits")
    cfgs = []
    for f in sorted(glob.glob(pat)):
        sdir = os.path.dirname(f)
        night = f.split(os.sep)[-4]
        cfg = os.path.basename(sdir).replace("keck_kcrm_", "")
        pf = glob.glob(os.path.join(sdir, "*.pypeit"))
        dk = "?"
        if pf:
            m = re.search(r"decker:\s*(\S+)", open(pf[0]).read())
            if m:
                dk = m.group(1)
        d = load(f)
        if d:
            d.update(tag=f"{night} {cfg}", decker=dk)
            cfgs.append(d)
    if not cfgs:
        sys.exit(f"no sens files matched {pat}")

    order = sorted({c["decker"] for c in cfgs})
    seen = {k: 0 for k in order}
    for c in sorted(cfgs, key=lambda c: (c["decker"], c["tag"])):
        c["color"] = SLICER_COLOR.get(c["decker"], "#4a3aa7")
        c["ls"] = STYLES[seen[c["decker"]] % len(STYLES)]
        seen[c["decker"]] += 1

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 8.5), sharex=True,
        gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08))
    fig.patch.set_facecolor("white")

    for ax in (ax1, ax2):
        ax.set_facecolor("white")
        ax.grid(True, color=INK3, alpha=0.22, linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(INK3)
        ax.tick_params(colors=INK2, labelsize=10)

    for c in cfgs:
        ax1.plot(c["w"], c["zp"], ".", color=c["color"], ms=1.4, alpha=0.16,
                 rasterized=True)
    for c in cfgs:
        ax1.plot(c["w"], c["zf"], c["ls"], color=c["color"], lw=2.0, alpha=0.95)
        ax2.plot(c["tw"], c["th"], c["ls"], color=c["color"], lw=2.0, alpha=0.95)

    # The band where every config is >=150 A clear of its own fit boundary --
    # the only wavelengths where all of them can be compared without edge bias.
    blo = max(c["lo"] for c in cfgs) + 150
    bhi = min(c["hi"] for c in cfgs) - 150
    if bhi > blo:
        for ax in (ax1, ax2):
            ax.axvspan(blo, bhi, color=INK3, alpha=0.13, lw=0, zorder=0)
        ax1.annotate(f"comparison band\n{blo:.0f}-{bhi:.0f} A\n(150 A clear of every fit edge)",
                     xy=(0.5 * (blo + bhi), 0.97), xycoords=("data", "axes fraction"),
                     ha="center", va="top", fontsize=8.5, color=INK2)

    ax1.set_ylabel("zeropoint (AB mag)", color=INK, fontsize=11)
    ax2.set_ylabel("throughput (%)", color=INK, fontsize=11)
    ax2.set_xlabel("wavelength (A)", color=INK, fontsize=11)
    ax2.set_ylim(bottom=0)

    # Group medians in the comparison band -- the finding, stated on the plot.
    if bhi > blo:
        for dk in order:
            v = []
            for c in cfgs:
                if c["decker"] != dk:
                    continue
                m = (c["tw"] >= blo) & (c["tw"] <= bhi)
                if m.any():
                    v.append(np.median(c["th"][m]))
            if v:
                ax2.annotate(f"{dk}  {np.mean(v):.1f}%",
                             xy=(bhi, np.mean(v)), xytext=(14, 0),
                             textcoords="offset points", va="center",
                             fontsize=10, color=INK,
                             bbox=dict(boxstyle="round,pad=0.28", fc="white",
                                       ec=SLICER_COLOR.get(dk, INK3), lw=1.4))

    handles = [Line2D([], [], color=c["color"], ls=c["ls"], lw=2.0,
                      label=f"{c['tag']}  ({c['decker']})") for c in
               sorted(cfgs, key=lambda c: (c["decker"], c["tag"]))]
    leg = ax1.legend(handles=handles, loc="lower right", ncol=2, fontsize=9,
                     frameon=True, framealpha=0.95, edgecolor=INK3,
                     labelcolor=INK2)
    leg.get_frame().set_facecolor("white")

    grating = os.path.basename(runroot.rstrip("/")).split()[0]
    ax1.set_title(f"{grating} sensitivity, all configs (IR, polyorder 15, "
                  f"full slice-coverage trim)\ncolour = slicer; drawn only "
                  f"inside each fit's good-pixel range",
                  color=INK, fontsize=12.5, pad=12)

    out = out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              f"{grating.lower()}_all_configs_sens_overlay_{suffix}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")
    for c in sorted(cfgs, key=lambda c: (c["decker"], c["tag"])):
        m = (c["tw"] >= blo) & (c["tw"] <= bhi) if bhi > blo else np.zeros_like(c["tw"], bool)
        v = np.median(c["th"][m]) if m.any() else np.nan
        print(f"   {c['tag']:16s} {c['decker']:7s} {c['lo']:.0f}-{c['hi']:.0f} A"
              f"   band median {v:5.1f}%")


if __name__ == "__main__":
    main()
