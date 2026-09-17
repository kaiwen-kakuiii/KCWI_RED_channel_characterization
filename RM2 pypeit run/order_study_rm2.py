#!/usr/bin/env python
"""Settle the zeropoint polynomial order on RM2, by RH3 section 6's test.

The usual metrics cannot choose a polynomial order.  Interior scatter falls
monotonically as order rises, because **that is also what a polynomial following
noise produces** -- measured on RM1, 0.042 -> 0.027 while night-to-night shape
RMS blew out from 0.7% to 4.9%.  The test that works is external: independent
nights of the same star through the same grating at the same central wavelength
must produce the same *shape*.  Transparency scales a curve up or down; it does
not move the peak.

RM2 can run that test on **three groups**, one of them a triple:

    cenwave 8950   2024-06-11 B, 2025-01-01 B, 2025-01-02 B     feige34, all three
    cenwave 8900   2024-01-04 B, 2024-04-30 B                   feige34, both
    cenwave 8850   2024-06-10 A, 2024-12-24 A                   feige110 vs feige 34

A group of three gives three pairwise comparisons, and the WORST of them is the
number to read: an order that holds two nights together and loses the third has
not held.

**8850 is a weaker test than the other two.**  Its two nights carry different
standard stars, so any error in either archival spectrum enters the comparison as
shape disagreement that no polynomial order can fix.  Read it as corroboration,
never as the deciding group.

## Why the inherited order is not safe to assume

Degrees of freedom per Angstrom is the quantity that transfers between gratings,
not the order itself:

    RL     order 15 over ~3400 A   =  227 A per DOF      tuned, worked
    RH3    order 15 over ~750 A    =   47 A per DOF      WRONG: peaks 180 A apart
    RH3    order  5 over ~750 A    =  125 A per DOF      chosen; peaks 35 A apart
    RM1    order  5 over ~1300 A   =  215 A per DOF      chosen on this test
    RM2    order  ? over ~1990 A

At RM2's ~1990 A fit range, order 5 is 332 A/DOF -- STIFFER than anything that
has worked so far, which argues for going up, to 7 (249) or 9 (199).  That is an
argument, not a measurement, and on RM1 the same arithmetic pointed at 7 while
the night-pair test said 5.  This script is what decides.

## What is measured

Per order, per group:

    shape RMS        scatter between two nights after dividing out a single
                     per-night grey factor, over the wavelengths both cover with
                     >150 A of margin from either fit edge.  Grey division is the
                     point: it removes transparency, which is real and is what
                     averaging nights is for, and leaves only shape, which is
                     the fitter's business
    peak separation  distance between the two nights' throughput maxima
    interior scatter mean over the pair of |measured - fitted| zeropoint, >25%
                     in from each end.  Reported to show it falling while shape
                     agreement degrades -- the trap, not the criterion
    A per DOF        fit range divided by (order + 1)

Usage:
    python order_study_rm2.py --orders 3 5 7 9 15        # run + analyse
    python order_study_rm2.py --analyse-only             # re-read what exists
"""

import argparse
import glob
import os
import re
import subprocess

import numpy as np
from astropy.io import fits

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))

# The groups the test turns on: same grating, same central wavelength, different
# night.  Anything not in a group cannot contribute a shape comparison, which is
# why 7750 and 9850 do not appear -- each is a single visit.
#
# 8850's two nights carry DIFFERENT STARS (feige110 and feige 34).  The
# comparison is still worth making, but an archival-spectrum error enters it as
# shape disagreement, so it corroborates and does not decide.
PAIRS = {8850: ["2024-06-10/keck_kcrm_A", "2024-12-24/keck_kcrm_A"],
         8900: ["2024-01-04/keck_kcrm_B", "2024-04-30/keck_kcrm_B"],
         8950: ["2024-06-11/keck_kcrm_B", "2025-01-01/keck_kcrm_B",
                "2025-01-02/keck_kcrm_B"]}
CROSS_STAR = {8850}
EDGE = 150.0


def sens_of(sdir, order):
    return sorted(glob.glob(os.path.join(sdir, "order_study",
                                         f"*_sens_IR_p{order}cov.fits")))


def curve(path):
    """(wave, throughput, lo, hi, interior scatter) inside the fitted range only.

    Outside SENS_ZEROPOINT_FIT_GPM the polynomial extrapolates and has been seen
    to exceed 100%.  THROUGHPUT is an ImageHDU on the WAVE grid, not a column of
    SENS and not on the SENS_WAVE grid -- pairing it with SENS_WAVE silently
    drops it (RH2.md 4.3).
    """
    with fits.open(path) as hdu:
        d = hdu["SENS"].data
        sw = np.asarray(d["SENS_WAVE"]).ravel()
        zp = np.asarray(d["SENS_ZEROPOINT"]).ravel()
        zf = np.asarray(d["SENS_ZEROPOINT_FIT"]).ravel()
        gpm = np.asarray(d["SENS_ZEROPOINT_FIT_GPM"]).ravel().astype(bool)
        names = [h.name for h in hdu]
        if "THROUGHPUT" not in names or "WAVE" not in names:
            return None
        thr = np.asarray(hdu["THROUGHPUT"].data).ravel()
        tw = np.asarray(hdu["WAVE"].data).ravel()

    ok = gpm & np.isfinite(sw) & (sw > 0) & np.isfinite(zp) & np.isfinite(zf)
    if not ok.any():
        return None
    lo, hi = sw[ok].min(), sw[ok].max()
    q = 0.25 * (hi - lo)
    inner = ok & (sw > lo + q) & (sw < hi - q)
    scat = np.std((zp - zf)[inner]) if inner.sum() > 5 else np.nan

    m = (tw >= lo) & (tw <= hi) & np.isfinite(thr) & (thr > 0)
    return dict(w=tw[m], t=thr[m] * 100.0, lo=lo, hi=hi, scatter=scat)


def compare(a, b):
    """Shape RMS and peak separation for two nights, on their common interior.

    A single grey factor per pair is divided out first -- the ratio of the two
    curves' means over the common window.  That is exactly the transparency term
    which is real and which averaging nights exists to handle; what survives is
    shape, which is the fitter's.
    """
    lo = max(a["lo"], b["lo"]) + EDGE
    hi = min(a["hi"], b["hi"]) - EDGE
    if hi - lo < 100:
        return None
    grid = np.linspace(lo, hi, 400)
    ta = np.interp(grid, a["w"], a["t"])
    tb = np.interp(grid, b["w"], b["t"])
    grey = ta.mean() / tb.mean()
    resid = ta - tb * grey
    return dict(rms=float(np.sqrt(np.mean(resid ** 2)) / ta.mean() * 100.0),
                peak_a=float(grid[np.argmax(ta)]), peak_b=float(grid[np.argmax(tb)]),
                dpeak=float(abs(grid[np.argmax(ta)] - grid[np.argmax(tb)])),
                grey=float(grey), lo=lo, hi=hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", type=int, nargs="+", default=[3, 5, 7, 9, 15])
    ap.add_argument("--analyse-only", action="store_true")
    args = ap.parse_args()

    if not args.analyse_only:
        for cw, cfgs in PAIRS.items():
            for cfg in cfgs:
                night, setup = cfg.split("/")
                sdir = os.path.join(BASE, night, "pypeit_run", setup)
                spec1d = sorted(glob.glob(os.path.join(sdir, "spec1d_*.fits")))
                base = sorted(glob.glob(os.path.join(sdir, "*_IR_p*cov.sens")))
                if not spec1d or not base:
                    print(f"{cfg}: no spec1d / .sens yet -- run "
                          f"run_rm2_throughput.py first, skipped", flush=True)
                    continue
                # Reuse the trim the main run measured, so ONLY the order varies.
                m = re.search(r"trim_std_pixs\s*=\s*(\d+),\s*(\d+)",
                              open(base[-1]).read())
                if not m:
                    print(f"{cfg}: no trim in {base[-1]}, skipped", flush=True)
                    continue
                blue, red = m.group(1), m.group(2)
                tag = os.path.basename(spec1d[-1])[len("spec1d_"):].split(".fits")[0]
                out = os.path.join(sdir, "order_study")
                os.makedirs(out, exist_ok=True)
                for o in args.orders:
                    sfile = os.path.join(out, f"{tag}_IR_p{o}cov.sens")
                    with open(sfile, "w") as fh:
                        fh.write("[sensfunc]\n    algorithm = IR\n    extr = BOX\n"
                                 f"    trim_std_pixs = {blue}, {red}\n"
                                 f"    polyorder = {o}\n")
                    dest = os.path.join(out, f"{tag}_sens_IR_p{o}cov.fits")
                    if os.path.exists(dest):
                        print(f"{cfg} order {o}: exists, skipped", flush=True)
                        continue
                    print(f"{cfg} order {o}: running sensfunc ...", flush=True)
                    # Serial on purpose: parallel sensfuncs race on the shared
                    # telluric cache in ~/.cache/pypeit.
                    rc = subprocess.run(
                        [os.path.join(ENV, "pypeit_sensfunc"),
                         os.path.basename(spec1d[-1]), "-s", sfile, "-o", dest],
                        cwd=sdir, capture_output=True, text=True).returncode
                    if rc:
                        print(f"   FAILED rc={rc}", flush=True)

    print("\n" + "=" * 78)
    print("shape agreement between two independent nights, per polynomial order")
    print("=" * 78)
    for cw, cfgs in PAIRS.items():
        note = ("   (different stars -- corroboration only)"
                if cw in CROSS_STAR else "")
        print(f"\ncenwave {cw}:  {len(cfgs)} nights{note}")
        # Every pairing inside the group, because an order that holds two nights
        # and loses the third has not held.  With two nights this is the single
        # comparison RM1 made; with three it is three, and the WORST is the one
        # that decides.
        pairs = [(a, b) for i, a in enumerate(cfgs) for b in cfgs[i + 1:]]
        print(f"  {'order':>6}{'A/DOF':>8}{'worstRMS':>10}{'worst d(pk)':>12}"
              f"{'interior':>10}   per-pair shapeRMS / d(peak)")
        for o in args.orders:
            curves = {}
            for cfg in cfgs:
                night, setup = cfg.split("/")
                sdir = os.path.join(BASE, night, "pypeit_run", setup)
                f = sens_of(sdir, o)
                curves[cfg] = curve(f[0]) if f else None
            if any(c is None for c in curves.values()):
                miss = [c for c in cfgs if curves[c] is None]
                print(f"  {o:>6}   (missing: {', '.join(miss)})")
                continue
            res = []
            for a, b in pairs:
                r = compare(curves[a], curves[b])
                if r is not None:
                    res.append((a, b, r))
            if not res:
                print(f"  {o:>6}   (no common interior)")
                continue
            worst_rms = max(r["rms"] for _, _, r in res)
            worst_dpk = max(r["dpeak"] for _, _, r in res)
            adof = np.mean([c["hi"] - c["lo"] for c in curves.values()]) / (o + 1)
            sc = np.nanmean([c["scatter"] for c in curves.values()])
            detail = "  ".join(
                f"{a.split('/')[0][5:]}v{b.split('/')[0][5:]} "
                f"{r['rms']:.2f}%/{r['dpeak']:.0f}A" for a, b, r in res)
            print(f"  {o:>6}{adof:>8.0f}{worst_rms:>9.2f}%{worst_dpk:>11.0f}A"
                  f"{sc:>10.4f}   {detail}")
    print("\nRead shapeRMS and d(peak), not interior scatter: a polynomial")
    print("following noise lowers interior scatter while moving the peak.")


if __name__ == "__main__":
    main()
