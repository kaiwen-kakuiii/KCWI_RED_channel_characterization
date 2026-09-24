#!/usr/bin/env python
"""Settle the zeropoint polynomial order on RM1, by RH3 section 6's test.

RH3 established that the usual metrics cannot choose a polynomial order.
Interior scatter falls monotonically as order rises, because **that is also what
a polynomial following noise produces**.  The test that works is external: two
independent nights of the same star through the same grating at the same central
wavelength must produce the same *shape*.  Night-to-night transparency scales a
curve up or down; it does not move the peak.

RM1 can run that test **twice**, which no earlier grating could:

    cenwave 6630   2023-12-10 C  and  2023-12-11 C
    cenwave 7510   2023-12-10 B  and  2023-12-11 B

Two pairs at two different grating angles.  If both pairs pick the same order,
the answer is a property of the fit and not of one accident.

## Why the inherited order is not safe to assume

Degrees of freedom per Angstrom is the quantity that transfers between gratings,
not the order itself:

    RL     order 15 over ~3400 A   =  227 A per DOF      tuned, worked
    RH2    order 15 over ~700 A    =   47 A per DOF      kept; RH3 doubts it
    RH3    order 15 over ~750 A    =   47 A per DOF      WRONG: peaks 180 A apart
    RH3    order  5 over ~750 A    =  125 A per DOF      chosen; peaks 35 A apart
    RM1    order  ? over ~1400 A

At RM1's ~1400 A fit range, order 7 is 175 A/DOF and order 5 is 233 A/DOF -- both
in RL's neighbourhood, and both far stiffer than what broke RH3.  That makes 7 a
reasonable default and not a settled answer, which is what this script is for.

## What is measured

Per order, per cenwave pair:

    shape RMS        scatter between the two nights after dividing out a single
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
    python order_study_rm1.py --orders 3 5 7 9 15        # run + analyse
    python order_study_rm1.py --analyse-only             # re-read what exists
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

# The pairs the test turns on: same star, same grating, same central wavelength,
# different night.  Anything not in a pair cannot contribute a shape comparison.
PAIRS = {6630: ["2023-12-10/keck_kcrm_C", "2023-12-11/keck_kcrm_C"],
         7510: ["2023-12-10/keck_kcrm_B", "2023-12-11/keck_kcrm_B"]}
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
                          f"run_rm1_throughput.py first, skipped", flush=True)
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
        print(f"\ncenwave {cw}:  {cfgs[0]}  vs  {cfgs[1]}")
        print(f"  {'order':>6}{'A/DOF':>8}{'shapeRMS':>10}{'d(peak)':>9}"
              f"{'peak A':>9}{'peak B':>9}{'grey':>7}{'interior':>10}")
        for o in args.orders:
            cs = []
            for cfg in cfgs:
                night, setup = cfg.split("/")
                sdir = os.path.join(BASE, night, "pypeit_run", setup)
                f = sens_of(sdir, o)
                cs.append(curve(f[0]) if f else None)
            if any(c is None for c in cs):
                print(f"  {o:>6}   (missing)")
                continue
            r = compare(cs[0], cs[1])
            if r is None:
                print(f"  {o:>6}   (no common interior)")
                continue
            adof = np.mean([c["hi"] - c["lo"] for c in cs]) / (o + 1)
            sc = np.nanmean([c["scatter"] for c in cs])
            print(f"  {o:>6}{adof:>8.0f}{r['rms']:>9.2f}%{r['dpeak']:>8.0f}A"
                  f"{r['peak_a']:>9.0f}{r['peak_b']:>9.0f}{r['grey']:>7.3f}"
                  f"{sc:>10.4f}")
    print("\nRead shapeRMS and d(peak), not interior scatter: a polynomial")
    print("following noise lowers interior scatter while moving the peak.")


if __name__ == "__main__":
    main()
