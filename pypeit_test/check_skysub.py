#!/usr/bin/env python
"""Check pypeit sky subtraction quality in a SlicerIFU spec2d file.

Usage:
    python check_skysub.py spec2d_*.fits [spec2d_*.fits ...]

WHAT IS MEASURED
----------------
Each IFU slice is split by the sky mask (user_regions, read from the .par file
next to Science/):

    slice width  0% ........................................... 100%
    |------------------|=====================|------------------|
    0%                 A                     B                100%
    <--- SKY ZONE ----->   <-- MASKED ZONE -->  <--- SKY ZONE --->

  SKY ZONE     pixels the sky model was FITTED on
  MASKED ZONE  pixels excluded from the fit (this is where the object sits)

The sky model is fitted on the sky zones but evaluated everywhere, so
"residual" (data - sky model) means different things in the two zones:

  residual in SKY ZONE     should be ~0. The fit was trying to match these
                           pixels, so a large value means the sky model is wrong.
  residual in MASKED ZONE  expected to be large -- that is the object, which
                           has no model (OBJMODEL=0 for IFU) and is measured
                           later at the datacube stage. Not an error.

VERDICT is therefore based on the SKY ZONE residual only.
"""
import sys
import os
import glob
import re
import numpy as np
from astropy.io import fits


def read_user_regions(spec2d_path):
    """Find 'user_regions = :A,B:' in the .par file beside Science/. Returns (A,B) or None."""
    rundir = os.path.dirname(os.path.dirname(os.path.abspath(spec2d_path)))
    pars = sorted(glob.glob(os.path.join(rundir, '*.par')), key=os.path.getmtime, reverse=True)
    for par in pars:                       # newest first; stale copies may say None
        m = re.search(r'user_regions\s*=\s*:\s*(\d+)\s*,\s*(\d+)\s*:', open(par).read())
        if m:
            return float(m.group(1)), float(m.group(2))
    return None


def report(path):
    h = fits.open(path, memmap=True)
    sci, sky, bpm = h[1].data, h[3].data, h[9].data
    med, std, slits = h[13].data, h[14].data, h[10].data
    good = bpm == 0
    nspat = sci.shape[1]
    xx = np.arange(nspat)[None, :]
    res = sci - sky

    reg = read_user_regions(path)
    print(f"\n=== {os.path.basename(path)} ===")
    if reg is None:
        print("  no user_regions found in .par -- cannot separate fitted vs masked zones.")
        print("  Falling back to whole-slice statistics (MED_CHIS / STD_CHIS only).")
    else:
        print(f"  sky mask: masked zone = {reg[0]:.0f}%-{reg[1]:.0f}% of each slice; "
              f"sky zones = 0-{reg[0]:.0f}% and {reg[1]:.0f}-100%")

    print(f"\n  {'slice':>5} {'STD_CHIS':>9} {'MED_CHIS':>9} "
          f"{'sky-zone resid':>15} {'masked-zone resid':>18}   verdict")
    bad_sky, leaks = [], []
    for i in range(len(std)):
        lo = int(np.median(slits['left_init'][i]))
        hi = int(np.median(slits['right_init'][i]))
        n = hi - lo
        inslit = (xx >= lo) & (xx <= hi) & good
        if inslit.sum() < 500:
            continue
        if reg is None:
            skyz = inslit
            maskz = None
        else:
            a, b = lo + int(reg[0] / 100 * n), lo + int(reg[1] / 100 * n)
            skyz = (((xx >= lo) & (xx < a)) | ((xx > b) & (xx <= hi))) & good
            maskz = ((xx >= a) & (xx <= b)) & good
        rs = np.median(res[skyz]) if skyz.sum() else np.nan
        rm = np.median(res[maskz]) if (maskz is not None and maskz.sum()) else np.nan
        # sky level for scale, from the sky zone
        lvl = max(np.median(sky[skyz]), 1.0)
        frac = abs(rs) / lvl
        if frac > 0.05:
            v = "SKY MODEL SUSPECT"
            bad_sky.append(i)
        else:
            v = "sky ok"
        if reg is not None and np.isfinite(rm) and rm > 5 * max(abs(rs), 1.0):
            leaks.append(i)
        print(f"  {i:5d} {std[i]:9.3f} {med[i]:9.3f} {rs:15.2f} "
              f"{rm if np.isfinite(rm) else float('nan'):18.2f}   {v}")

    print()
    if bad_sky:
        print(f"  WARNING: sky model suspect in {len(bad_sky)} slice(s): {bad_sky}")
        print(f"           sky-zone residual exceeds 5% of the sky level there, i.e. the model")
        print(f"           does not reproduce the pixels it was fitted on.")
    else:
        print("  Sky model reproduces the fitted sky zones in every slice (residual < 5% of sky).")
    if leaks:
        print(f"  Note: slices {leaks} show flux in the masked zone well above the sky-zone")
        print( "        residual. Expected where the object sits; if it extends to the mask")
        print( "        edges some object flux may lie outside the mask.")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        report(p)
