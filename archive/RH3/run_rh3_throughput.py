#!/usr/bin/env python
"""Take the RH3 spec2d files through to throughput, per THROUGHPUT_PROCEDURE.md.

    spec2d  --pypeit_coadd_datacube-->  datacube + whitelight
            --pypeit_extract_datacube-> spec1d (BOX aperture; -b ARCSEC)
            --measure_trim.py --spec2d-> trim_std_pixs
            --pypeit_sensfunc---------> zeropoint + throughput

## What RH3 is

Four nights, one usable configuration each (setup B; setup A is the night's stray
DARK plus biases every time), **Medium slicer, 2x2 binning** throughout:

    2023-11-12 B   cenwave 8600   feige110   2 frames  (180, 360 s)
    2024-10-01 B   cenwave 8400   g191b2b    3 frames  (100 s each; a 1 s
                                             acquisition frame is dropped)
    2024-10-29 B   cenwave 8600   feige110   2 frames  (60, 300 s)
    2024-11-06 B   cenwave 8600   feige110   2 frames  (300 s each)

Unlike RH4 this is **one star per configuration**, so the TARGNAME split that
run_rh4_throughput.py exists for never fires here.  It is kept rather than
removed: it costs nothing, and a driver that silently coadds two stars is the
failure RH4 documents.  What RH3 gains instead is the cross-check RH2 got --
g191b2b on 2024-10-01 against feige110 on the other three nights, two different
standards through the same grating.

## Two central wavelengths, and they are not interchangeable

Three nights sit at cenwave 8600 and one at 8400.  RH2 established that a
high-dispersion KCRM grating's throughput peak tracks the grating angle (mean
|peak - cenwave| = 38 A for RH2 against 1343 A for RL), so the 8400 config
samples a different point of the blaze and must be read as its own curve, not
averaged with the 8600 ones at a fixed wavelength.  `plot_rh1_throughput.py -g
RH3` colours by cenwave for exactly this reason.

## One bad pixel can truncate a night, so the trim bridges short gaps

`measure_trim.py` keeps the largest CONTIGUOUS block of pixels above the counts
threshold, which exists to stop an isolated cosmic ray anchoring the range far
out into dead spectrum.  RH3 hit the same mechanism running the other way: on
2024-10-29 one pixel at 8661.8 A reading -184.7 counts split an otherwise
continuous spectrum into a 1296-px blue block and an 872-px red block, and the
blue one won:

    bridge 0   keep 8208.5-8657.3 A   trim 176, 935     <- 353 A thrown away
    bridge 5   keep 8208.5-8964.5 A   trim 176, 176

The discarded region carries 270-520 counts against a 107-count threshold, so it
is good data.  `--bridge 5` joins runs separated by up to five pixels first.  It
changes nothing on the other three RH3 nights, and `measure_trim`'s own default
stays 0 so no earlier grating's numbers move.

## Aperture

`-b` is set in ARCSEC and fixed across cubes, for the reason RH1 spent a day on
(RH1_PROCEDURE.md problem 5): `pypeit_extract_datacube`'s default aperture is
4*sigma with sigma fit in SPAXELS, and a spaxel is 1.358" on Large/2x2 but
0.679" on Medium/2x2, so the default measures a different solid angle per
configuration.  **RH3 is Medium/2x2, not Large like RH1 and RH4** -- which is
precisely why the radius has to be stated in arcsec: 3.4" is the same solid
angle here as there, while the default would silently halve it.

Usage:
    python run_rh3_throughput.py [--boxcar 3.4] [--max-dam 0.10] [--polyorder 15]
    python run_rh3_throughput.py --boxcar 3.4 --skip-coadd     # re-extract only
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
from astropy.io import fits

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
TRIM = os.path.join(ROOT, "pypeit_test", "measure_trim.py")


def groups(spec2ds, max_dam, gap_min=None):
    """Split frames into cubes: by STAR first, then by airmass and pointing.

    Star is a hard split and comes first -- see the module docstring.  Within one
    star the RH1 criteria apply unchanged: airmass is what sets the extinction
    correction, and pypeit applies that correction once per cube from a single
    header airmass (coadd3d.py disables the per-frame correction outright).
    Measured against PypeIt's Mauna Kea curve, an 0.10 airmass spread is ~0.8% in
    flux, under the ~1% scatter this measurement resolves.
    """
    rows = []
    for f in spec2ds:
        h = fits.getheader(f)
        star = str(h.get("TARGNAME") or "").strip().lower()
        star = re.sub(r"[^a-z0-9]+", "", star) or "std"
        rows.append((star, h.get("MJD", 0.0), h.get("AIRMASS", 0.0),
                     (round(h.get("RAOFF", 0.0) or 0.0, 2),
                      round(h.get("DECOFF", 0.0) or 0.0, 2)), f))
    rows.sort()

    out = []
    for star in sorted({r[0] for r in rows}):
        sub = sorted([r for r in rows if r[0] == star], key=lambda r: r[1])
        cur = [sub[0]]
        for prev, this in zip(sub, sub[1:]):
            ams = [r[2] for r in cur] + [this[2]]
            split = (this[3] != cur[0][3]                    # pointing changed
                     or max(ams) - min(ams) > max_dam)       # airmass drifted
            if gap_min is not None and (this[1] - prev[1]) * 1440.0 > gap_min:
                split = True
            if split:
                out.append((star, cur))
                cur = []
            cur.append(this)
        out.append((star, cur))
    return out


def sh(cmd, log, cwd):
    with open(log, "a") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def check_cube(cubepath, tag, boxcar=None, spaxel=None):
    """Cube geometry guard, carried over from post_rh2.py.

    Rejects only a degree-scale CDELT -- the documented single-frame failure that
    inflates extracted counts by ~1e21 and ruins the sensfunc silently.  An
    oversized spatial grid is NOT rejected: three RH2 configs came out 553-1307
    spaxels on one axis and were measured to be fine (99.82% of the flux inside
    16 rows), and an earlier version that rejected on field size would have
    thrown away every valid RH2 cube.  Poor flux concentration only warns.

    Concentration is measured inside the EXTRACTION APERTURE, not inside a fixed
    +/-8 spaxels.  RH2's grids were hundreds of spaxels across, so +/-8 was a
    core; an RH4 cube is 15 x 23, where +/-8 spans nearly the whole field and the
    metric ends up measuring the cube's edge ring instead of the star.  Measured
    on RH4's feige110_2024-12-28_B: 77.6% within +/-8 spaxels, reading as a smeared
    stack, while the star is in fact textbook -- peak at spaxel (7,11), the grid
    centre and the same spaxel as the known-good RL 2024-12-04 B cube, with 63%
    of the signed total inside ONE spaxel.  The missing flux is a spatial
    resampling artifact confined to the outermost row and column (rows 0/14 hold
    10.7%/6.9%, columns 0/22 hold 11.5%/9.8%), 4-11 spaxels from the peak and so
    entirely outside a 3.4" aperture.  The edge ring is reported separately
    rather than folded into the concentration number.
    """
    try:
        with fits.open(cubepath) as hdu:
            hdr = hdu[1].header if len(hdu) > 1 else hdu[0].header
            cd3 = np.asarray(hdu[1].data, float)
        shape = cd3.shape
        cd = [abs(float(hdr.get(f"CDELT{i}", np.nan))) for i in (1, 2)]
        nx = [int(hdr.get(f"NAXIS{i}", 0)) for i in (1, 2)]
        ext = [n * c * 3600.0 for n, c in zip(nx, cd) if np.isfinite(c)]
        bad_wcs = any(c > 0.01 for c in cd if np.isfinite(c))

        core, ring, rad = np.nan, np.nan, np.nan
        wl = np.nansum(np.where(np.isfinite(cd3), cd3, 0), axis=0)
        tot = wl.sum()
        if tot > 0:
            j = np.unravel_index(np.argmax(wl), wl.shape)
            # Aperture radius in spaxels, from the arcsec aperture actually used.
            sp = spaxel if spaxel else (cd[1] * 3600.0 if np.isfinite(cd[1]) else 1.0)
            rad = (boxcar / sp) if (boxcar and sp > 0) else 2.5
            yy, xx = np.mgrid[0:wl.shape[0], 0:wl.shape[1]]
            r = np.hypot(yy - j[0], xx - j[1])
            core = wl[r <= rad].sum() / tot
            edge = np.zeros(wl.shape, bool)
            edge[0, :] = edge[-1, :] = True
            edge[:, 0] = edge[:, -1] = True
            ring = wl[edge].sum() / tot
        del cd3, wl

        print(f"  cube {shape}  CDELT1/2 = {cd[0]:.3g}/{cd[1]:.3g} deg  "
              f"field {ext[0]:.0f}\" x {ext[1]:.0f}\"  "
              f"flux within {rad:.1f} spaxels of peak {100*core:.1f}%  "
              f"(outer ring {100*ring:.1f}%)"
              + ("   *** WCS SCALE BROKEN ***" if bad_wcs else ""), flush=True)
        if bad_wcs:
            print("  refusing to extract: a degree-scale CDELT makes counts "
                  "meaningless", flush=True)
            return False
        if np.isfinite(core) and core < 0.60:
            print(f"  WARNING star is not compact ({100*core:.1f}% inside the "
                  f"aperture) -- check the whitelight for a smeared or doubled "
                  f"stack", flush=True)
        if np.isfinite(ring) and ring > 0.15:
            print(f"  NOTE {100*ring:.1f}% of the signed flux sits in the cube's "
                  f"outermost row/column -- a spatial resampling edge artifact, "
                  f"outside the aperture and harmless to the extraction",
                  flush=True)
    except Exception as e:
        print(f"  could not read cube ({type(e).__name__}), continuing", flush=True)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-dam", type=float, default=0.10,
                    help="largest airmass spread allowed inside one cube")
    ap.add_argument("--gap", type=float, default=None,
                    help="optional extra split on a time gap, in minutes")
    ap.add_argument("--polyorder", type=int, default=15)
    ap.add_argument("--boxcar", type=float, default=3.4,
                    help="extraction radius in ARCSEC.  Fixed across cubes on "
                         "purpose; the default aperture is 4*sigma with sigma "
                         "fit in spaxels, which is a different solid angle per "
                         "configuration")
    ap.add_argument("--bridge", type=int, default=5,
                    help="join above-threshold runs separated by at most this "
                         "many pixels before taking the largest.  RH3 needs it: "
                         "on 2024-10-29 a SINGLE pixel at 8661.8 A reading -184.7 "
                         "counts split the spectrum in two, and the blue half won "
                         "on length, discarding 353 A that carries 270-520 counts "
                         "against a 107-count threshold.  Verified to change "
                         "nothing on the other RH3 nights; measure_trim's own "
                         "default is 0, so earlier gratings are untouched")
    ap.add_argument("--only", default=None)
    ap.add_argument("--skip-coadd", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for sdir in sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run", "keck_kcrm_*"))):
        if args.only and args.only not in sdir:
            continue
        s2d = sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        if not s2d:
            continue
        night = os.path.relpath(sdir, BASE).split(os.sep)[0]
        letter = os.path.basename(sdir).replace("keck_kcrm_", "")
        log = os.path.join(sdir, "throughput.log")

        grps = groups(s2d, args.max_dam, gap_min=args.gap)
        percount = {}
        for star, _ in grps:
            percount[star] = percount.get(star, 0) + 1

        seen = {}
        for star, grp in grps:
            seen[star] = seen.get(star, 0) + 1
            am = np.array([r[2] for r in grp])
            mjd = np.array([r[1] for r in grp])
            tag = f"{star}_{night}_{letter}"
            if percount[star] > 1:
                tag += f"_v{seen[star]}"
            print(f"\n=== {tag}: {len(grp)} frames, airmass {am.min():.2f}-"
                  f"{am.max():.2f} (spread {am.max()-am.min():.2f}), "
                  f"span {(mjd.max()-mjd.min())*1440:.0f} min", flush=True)

            cube = f"{tag}.fits"
            c3d = os.path.join(sdir, f"{tag}.coadd3d")
            with open(c3d, "w") as fh:
                fh.write(f"# {star}, {night} setup {letter}; airmass "
                         f"{am.min():.2f}-{am.max():.2f}\n"
                         "# align is deliberately NOT set: measured on RH2 "
                         "2023-11-07 B, align=True and align=False give\n"
                         "# byte-identical cube geometry, because the header WCS "
                         "already places a small nod correctly.\n"
                         "[rdx]\n    spectrograph = keck_kcrm\n\n[reduce]\n"
                         "    [[cube]]\n        combine = True\n"
                         f"        output_filename = {cube}\n"
                         "        save_whitelight = True\n\nspec2d read\nfilename\n")
                for r in grp:
                    fh.write(f"Science/{os.path.basename(r[-1])}\n")
                fh.write("spec2d end\n")
            if args.dry_run:
                for r in grp:
                    print(f"    {os.path.basename(r[-1])}  airmass {r[2]:.2f}")
                continue

            cubepath = os.path.join(sdir, cube)
            if args.skip_coadd and os.path.exists(cubepath):
                print("  reusing existing cube", flush=True)
            elif sh([os.path.join(ENV, "pypeit_coadd_datacube"),
                     os.path.basename(c3d), "-o"], log, sdir):
                print(f"  COADD FAILED -- see {log}", flush=True)
                continue
            if not os.path.exists(cubepath):
                print(f"  COADD produced no cube -- see {log}", flush=True)
                continue
            if not check_cube(cubepath, tag, boxcar=args.boxcar):
                continue

            xcmd = [os.path.join(ENV, "pypeit_extract_datacube"), cube, "-o"]
            if args.boxcar is not None:
                xcmd += ["-b", str(args.boxcar)]
            if sh(xcmd, log, sdir):
                print(f"  EXTRACT FAILED -- see {log}", flush=True)
                continue
            spec1d = sorted(glob.glob(os.path.join(sdir, f"spec1d_{tag}*.fits")))
            if not spec1d:
                print("  no spec1d produced", flush=True)
                continue

            # --spec2d restricts the fit to wavelengths ALL slits recorded.
            # Without it the cube's partial-coverage ends are included, where the
            # coadd fabricates flux by resampling across missing slices -- and
            # the counts look completely normal there.
            out = subprocess.run(
                [os.path.join(ENV, "python"), TRIM, "--spec2d",
                 os.path.join(sdir, "Science"),
                 "--bridge", str(args.bridge), spec1d[-1]],
                capture_output=True, text=True)
            open(log, "a").write(out.stdout + out.stderr)
            m = re.search(r"trim_std_pixs = (\d+),\s*(\d+)", out.stdout)
            if not m:
                print(f"  TRIM FAILED -- see {log}", flush=True)
                continue
            blue, red = m.group(1), m.group(2)

            sfile = os.path.join(sdir, f"{tag}_IR_p{args.polyorder}cov.sens")
            with open(sfile, "w") as fh:
                fh.write("# IR sensfunc; trim = counts criterion AND full\n"
                         "# slice-coverage interval (measure_trim.py --spec2d).\n"
                         "[sensfunc]\n    algorithm = IR\n    extr = BOX\n"
                         f"    trim_std_pixs = {blue}, {red}\n"
                         f"    polyorder = {args.polyorder}\n")
            sens = f"{tag}_sens_IR_p{args.polyorder}cov.fits"
            # Strictly one at a time: parallel sensfuncs race on the shared
            # telluric cache in ~/.cache/pypeit and one dies on a half-written
            # entry (THROUGHPUT_PROCEDURE.md).
            rc = sh([os.path.join(ENV, "pypeit_sensfunc"),
                     os.path.basename(spec1d[-1]), "-s", os.path.basename(sfile),
                     "-o", sens], log, sdir)
            print(f"  trim {blue},{red}  -> "
                  f"{'sensfunc OK' if rc == 0 else f'SENSFUNC FAILED rc={rc}'}: {sens}",
                  flush=True)
            if rc == 0:
                summarise(os.path.join(sdir, sens), tag)

    print("\nALL SETUPS PROCESSED", flush=True)


def summarise(path, tag):
    """Fit quality and throughput, read only inside the fitted range.

    Outside SENS_ZEROPOINT_FIT_GPM the polynomial extrapolates and can be wrong
    by magnitudes.  Interior/edge windows are PROPORTIONAL, not the fixed 300 A
    RL used: RH2 showed a fixed window leaves ~90 A on a high-dispersion config
    and reports a meaningless number.  THROUGHPUT is ImageHDU 5, not a column of
    the SENS table.
    """
    with fits.open(path) as hdu:
        d = hdu["SENS"].data
        w = np.asarray(d["SENS_WAVE"]).ravel()
        zp = np.asarray(d["SENS_ZEROPOINT"]).ravel()
        zf = np.asarray(d["SENS_ZEROPOINT_FIT"]).ravel()
        gpm = np.asarray(d["SENS_ZEROPOINT_FIT_GPM"]).ravel().astype(bool)
        # THROUGHPUT is an ImageHDU beside WAVE and ZEROPOINT, NOT a column of
        # the SENS table -- reading d["THROUGHPUT"] finds nothing and fails
        # silently (RH2.md section 4.3).  It also sits on the WAVE HDU's grid,
        # which is NOT the SENS_WAVE grid: measured here, 8408 points against
        # 10197.  Pairing THROUGHPUT with SENS_WAVE silently drops it.
        names = [h.name for h in hdu]
        thr = (np.asarray(hdu["THROUGHPUT"].data).ravel()
               if "THROUGHPUT" in names else None)
        twave = (np.asarray(hdu["WAVE"].data).ravel()
                 if "WAVE" in names else None)

    ok = gpm & np.isfinite(w) & (w > 0) & np.isfinite(zp) & np.isfinite(zf)
    if not ok.any():
        print(f"  {tag}: sens written but no good fit pixels", flush=True)
        return
    ww, resid = w[ok], zp[ok] - zf[ok]
    lo, hi = ww.min(), ww.max()
    q = 0.25 * (hi - lo)
    interior = (ww > lo + q) & (ww < hi - q)
    edge = (ww < lo + 0.05 * (hi - lo)) | (ww > hi - 0.05 * (hi - lo))
    sc = np.std(resid[interior]) if interior.sum() > 5 else np.nan
    ee = np.max(np.abs(resid[edge])) if edge.sum() else np.nan
    print(f"  fit range {lo:.0f}-{hi:.0f} A   interior scatter {sc:.4f} mag   "
          f"worst edge error {ee:.3f} mag", flush=True)
    if thr is not None and twave is not None and thr.size == twave.size:
        # Restrict to the fitted range by WAVELENGTH, since the two grids differ.
        m = (twave >= lo) & (twave <= hi) & np.isfinite(thr) & (thr > 0)
        if m.any():
            tw, tv = twave[m], thr[m]
            q = 0.25 * (tw.max() - tw.min())
            inner = (tw > tw.min() + q) & (tw < tw.max() - q)
            print(f"  throughput interior mean {100*np.mean(tv[inner]):.1f}%   "
                  + "  ".join(f"{p}A {100*tv[np.argmin(np.abs(tw-p))]:.1f}%"
                              for p in (8200, 8400, 8600, 8800, 9000) if lo < p < hi),
                  flush=True)


if __name__ == "__main__":
    main()
