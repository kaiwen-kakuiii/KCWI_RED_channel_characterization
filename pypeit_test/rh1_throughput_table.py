#!/usr/bin/env python
"""Collect the RH1 throughput results, and check each cube stacked properly.

Two things are reported per cube, because a sensible-looking throughput curve
can still come from a cube whose exposures failed to stack:

  whitelight star  the collapsed image should hold ONE compact source.  Smeared
                   or doubled means the frames did not align, and every count
                   the aperture sums afterwards is suspect (THROUGHPUT_PROCEDURE
                   step 1).  Reported as `encl`, the fraction of total flux the
                   extraction aperture actually enclosed.

                   Measured in ARCSEC, not spaxels.  The spaxel is 1.358" on
                   Large/2x2 but 0.339" on Small/1x1, so a radius quoted in
                   spaxels means a different solid angle per night and the
                   column would not be comparable -- the same units trap that
                   made the default 4-sigma aperture unusable across slicers.

  throughput       read ONLY inside SENS_ZEROPOINT_FIT_GPM.  Outside the fitted
                   range the polynomial extrapolates and has been seen to give
                   throughput above 100%.

Nights in DROP are excluded from every mean; see that constant for why.

Usage:
    python rh1_throughput_table.py [<run dir>]
"""
import collections
import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "RH1 pypeit run")
SAMPLE = [6200, 6400, 6600, 6800]
EDGE = 150.0   # A; within this of a fit boundary the zeropoint fit is unreliable
               # (THROUGHPUT_PROCEDURE: ~3% interior, up to ~12% near an edge)

# Nights dropped from the means because the photons did not arrive.  Compared
# within one star, so expected flux is identical, over 6050-6350 A with
# extinction removed:
#
#     g191b2b  2023-09-16   932 cts/s -> 30.2%      2023-11-07   793 -> 25.3%
#     feige34  2023-11-17  1660 cts/s -> 29.6%      2024-03-13  1170 -> 20.8%
#
# The deficit is already in the extracted counts, so sensfunc, the zeropoint
# polynomial and the trim interval did not cause it.  Aperture accounts for only
# part: corrected to total flux the two still sit 14% and 22% low.  What is left
# is transparency / focus / telescope state, which these files cannot separate.
# The nights that behave agree at ~30% on two different standards, matching the
# 29.6% found at 6400 and 6600 A -- so the blue end is not anomalous, the sample
# was.  Averaging these in pulled @6200 A down about 3 points.
DROP = {"2023-11-07": "non-photometric: 15% low in counts vs 2023-09-16, same star",
        "2024-03-13": "non-photometric: 30% low in counts vs 2023-11-17, same star"}


def pixscale(cube):
    """Cube spaxel size in arcsec, from whichever WCS keyword is present."""
    with fits.open(cube) as h:
        for hdu in h:
            for k in ("CD2_2", "CDELT2"):
                if k in hdu.header:
                    return abs(hdu.header[k]) * 3600.0
    return np.nan


def whitelight_quality(path, radius_arcsec, scale):
    """Fraction of flux the extraction aperture enclosed, and any rival peak.

    `radius_arcsec` is the aperture actually used, read back from the spec1d;
    `scale` converts it to spaxels for this particular cube's grid.
    """
    img = np.nan_to_num(fits.getdata(path))
    if img.ndim > 2:
        img = img[0]
    # Background-subtract and clip: the whitelight image carries slightly
    # negative sky residuals, and dividing by a total that includes them lets
    # the enclosed fraction exceed 1, which is meaningless.
    img = np.clip(img - np.median(img), 0, None)
    tot = img.sum()
    if tot <= 0 or not np.isfinite(scale) or scale <= 0:
        return np.nan, np.nan
    iy, ix = np.unravel_index(np.argmax(img), img.shape)
    yy, xx = np.mgrid[: img.shape[0], : img.shape[1]]
    r = np.hypot(yy - iy, xx - ix) * scale            # arcsec, not spaxels
    frac = float(img[r <= radius_arcsec].sum() / tot)
    far = img.copy()
    far[r <= 3 * scale] = 0
    sep = float(np.hypot(*(np.array(np.unravel_index(np.argmax(far), far.shape))
                           - np.array([iy, ix]))) * scale)
    return frac, sep


def main():
    rows = []
    for sens in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run",
                                              "keck_kcrm_*", "*_sens_*.fits"))):
        tag = os.path.basename(sens).split("_sens_")[0]
        sdir = os.path.dirname(sens)
        with fits.open(sens) as h:
            t = Table(h["SENS"].data)
            wave = np.asarray(h["WAVE"].data).flatten()
            thru = np.asarray(h["THROUGHPUT"].data).flatten()
        gpm = np.asarray(t["SENS_ZEROPOINT_FIT_GPM"][0], dtype=bool)
        sw = np.asarray(t["SENS_WAVE"][0]).flatten()
        lo, hi = sw[gpm].min(), sw[gpm].max()

        ok = (wave >= lo) & (wave <= hi) & np.isfinite(thru) & (thru > 0)
        vals, edge = {}, {}
        for lam in SAMPLE:
            inside = ok.any() and lo <= lam <= hi
            vals[lam] = float(np.interp(lam, wave[ok], thru[ok])) * 100 if inside else np.nan
            edge[lam] = inside and (lam - lo < EDGE or hi - lam < EDGE)
        mid = 0.5 * (lo + hi)
        centre = float(np.interp(mid, wave[ok], thru[ok])) * 100 if ok.any() else np.nan

        pf = glob.glob(os.path.join(sdir, "*.pypeit"))
        cw = np.nan
        if pf:
            mm = re.search(r"cenwave:\s*([0-9.]+)", open(pf[0]).read())
            if mm:
                # same nominal setting differs in the 2nd decimal between nights
                cw = round(float(mm.group(1)))
        cube = os.path.join(sdir, f"{tag}.fits")
        am = fits.getheader(cube).get("AIRMASS", np.nan) if os.path.exists(cube) else np.nan
        scale = pixscale(cube) if os.path.exists(cube) else np.nan
        s1 = os.path.join(sdir, f"spec1d_{tag}.fits")
        rad = (fits.getheader(s1, 1)["BOX_R_PIX"] * scale
               if os.path.exists(s1) and np.isfinite(scale) else np.nan)
        wl = glob.glob(os.path.join(sdir, f"{tag}_whitelight.fits"))
        frac, sep = (whitelight_quality(wl[0], rad, scale)
                     if wl and np.isfinite(rad) else (np.nan, np.nan))
        rows.append((tag, lo, hi, vals, frac, sep, edge, mid, centre, am, cw, rad))

    if not rows:
        sys.exit("no sensfunc files yet")

    def dropped(tag):
        return next((n for n in DROP if n in tag), None)

    print("values in (parentheses) sit within "
          f"{EDGE:.0f} A of a fit boundary -- unreliable, excluded from means")
    print("rows marked DROP are excluded from every mean -- see DROP in this file\n")
    print(f"{'cube':30s} {'cenw':>6s} {'airm':>5s} {'fit range (A)':>16s} " +
          " ".join(f"{l:>8d}" for l in SAMPLE) +
          f"{'centre':>16s} {'ap\"':>5s} {'encl':>5s}")
    for tag, lo, hi, vals, frac, sep, edge, mid, centre, am, cw, rad in rows:
        cells = []
        for l in SAMPLE:
            if not np.isfinite(vals[l]):
                cells.append("       -")
            elif edge[l]:
                cells.append(f"({vals[l]:5.1f}%)")
            else:
                cells.append(f" {vals[l]:5.1f}% ")
        c = f"{centre:5.1f}% @{mid:6.0f}" if np.isfinite(centre) else "        -"
        flag = ""
        if dropped(tag):
            flag = "  <-- DROP"
        elif np.isfinite(frac) and frac < 0.85:
            flag = "  <-- aperture missed flux: check stacking / focus"
        print(f"{tag:30s} {cw:6.0f} {am:5.2f} {lo:7.1f}-{hi:7.1f} " + " ".join(cells) +
              f" {c:>15s} {rad:5.2f} {frac:5.2f}{flag}")

    keep = [r for r in rows if not dropped(r[0])]
    print(f"\n  {len(rows) - len(keep)} of {len(rows)} cubes dropped: " +
          ", ".join(f"{n} ({why})" for n, why in DROP.items()))
    for lam in SAMPLE:
        v = np.array([r[3][lam] for r in keep if np.isfinite(r[3][lam]) and not r[6][lam]])
        if v.size:
            print(f"  @{lam} A   mean {v.mean():5.2f}%   std {v.std():4.2f}%   n={v.size} (interior only)")
    cen = np.array([r[8] for r in keep if np.isfinite(r[8])])
    if cen.size:
        print(f"  at each cube's own fit centre: mean {cen.mean():5.2f}%  "
              f"std {cen.std():4.2f}%  n={cen.size}")

    # Airmass can only be tested WITHIN a setup.  Across setups the central
    # wavelength differs, which rotates the grating and moves the blaze peak, so
    # even the same wavelength sits at a different point on the efficiency curve
    # -- and each cube's fit centre is a different wavelength again.  Comparing
    # visits of one setup holds cenwave, grating angle and star fixed and leaves
    # airmass as the only difference.
    byset = collections.defaultdict(list)
    for r in keep:
        base = re.sub(r"_v\d+$", "", r[0])
        byset[base].append(r)
    pairs = {k: v for k, v in byset.items() if len(v) > 1}
    if pairs:
        print("\nwithin-setup visit pairs (same cenwave, airmass is the only difference):")
        for k, v in sorted(pairs.items()):
            v = sorted(v, key=lambda r: r[9])
            lo_r, hi_r = v[0], v[-1]
            for lam in SAMPLE:
                a, b = lo_r[3][lam], hi_r[3][lam]
                if np.isfinite(a) and np.isfinite(b) and not lo_r[6][lam] and not hi_r[6][lam]:
                    print(f"  {k:28s} @{lam} A   airmass {lo_r[9]:.2f} -> {a:5.1f}%   "
                          f"airmass {hi_r[9]:.2f} -> {b:5.1f}%   "
                          f"diff {b-a:+5.1f} pts ({100*(b/a-1):+5.1f}%)")
    else:
        print("\nno completed visit pair yet -- the clean airmass test needs both "
              "visits of one setup")


if __name__ == "__main__":
    main()
