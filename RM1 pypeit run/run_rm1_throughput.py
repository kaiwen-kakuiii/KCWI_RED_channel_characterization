#!/usr/bin/env python
"""Take the RM1 spec2d files through to throughput, per THROUGHPUT_PROCEDURE.md.

    spec2d  --pypeit_coadd_datacube-->  datacube + whitelight
            --pypeit_extract_datacube-> spec1d (BOX aperture; -b ARCSEC)
            --measure_trim.py --spec2d-> trim_std_pixs
            --pypeit_sensfunc---------> zeropoint + throughput

## What RM1 is, and how it differs from every earlier run

Six nights, nine science configurations, **one star throughout** (feige34):

    2024-04-01 A   cenwave 6130   Small 1x1   2 frames  ( 25, 100 s)
    2024-05-02 B   cenwave 6200   Large 2x2   6 frames  (5, 15, 45, 60, 120, 120 s)
    2024-03-15 B   cenwave 6300   Large 2x2   7 frames  (6, 15, 45 x5)
    2023-12-09 B   cenwave 6480   Large 2x2   2 frames  (  7,  14 s)
    2023-12-10 C   cenwave 6630   Large 2x2   4 frames  (15 s each)
    2023-12-11 C   cenwave 6630   Large 2x2   4 frames  (15 s each)
    2024-04-01 B   cenwave 7390   Small 1x1   2 frames  (200 s each)
    2023-12-10 B   cenwave 7510   Large 2x2   3 frames  (15 s each)
    2023-12-11 B   cenwave 7510   Large 2x2   3 frames  (15 s each)

**Seven central wavelengths and two slicers.**  RL had four cenwaves, RH3 two,
RH2 three.  Nothing here can be read at "a" wavelength across configs without
saying which grating angle it came from -- RM1 spans 6130-7510 A of central
wavelength, 1380 A of grating rotation, the widest of any run in this project.

**Two configurations share a night** on 2023-12-10 and 2023-12-11 (setup B at
7510, setup C at 6630), so a night is not a configuration here and the driver
must key on the setup directory, never on the date.

**Two cenwaves have two independent nights each** -- 6630 and 7510, both on
2023-12-10 and 2023-12-11.  That is the pair RH3 section 6 needs to settle
polynomial order: three nights of the same star through the same grating at the
same central wavelength must produce the same SHAPE, and disagreement is the
fitter, not the sky.  RM1 gets that test twice over, at two different angles.
See `order_study_rm1.py`; do not accept the inherited polyorder without it.

## The aperture is 3.4 arcsec, and here that is not optional

`pypeit_extract_datacube`'s default aperture is 4*sigma with sigma fit in
SPAXELS, and a spaxel is **1.358" on Large/2x2 but 0.339" on Small/1x1** -- a
factor of four.  RM1 mixes both in one run (2024-04-01 is Small/1x1, the other
five nights Large/2x2), so the default would measure two different solid angles
and the Small night would read low, which is exactly the 2.4-point error
RH1_PROCEDURE.md problem 5 documents.  3.4" is the radius RH1 validated across
these same two configurations, where every cube came out 100% on-field.

## Trim bridging is kept on

`--bridge 5` joins above-threshold runs separated by at most five pixels before
the largest is taken.  RH3 needed it because ONE pixel reading -184.7 counts
split a spectrum in two and the blue half won on length, silently discarding
353 A of good data with nothing in the log to say so.  It changed nothing on the
nights that did not need it, so it stays on here; measure_trim's own default is
still 0 and no earlier grating moves.

Usage:
    python run_rm1_throughput.py [--boxcar 3.4] [--polyorder N] [--dry-run]
    python run_rm1_throughput.py --skip-coadd          # re-extract only
    python run_rm1_throughput.py --only 2024-03-15
"""
import argparse
import glob
import os
import re
import subprocess

import numpy as np
from astropy.io import fits
from scipy import ndimage

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
TRIM = os.path.join(ROOT, "pypeit_test", "measure_trim.py")

# Where to quote throughput in the per-config summary.  RM1's union runs
# ~5400-8230 A across the seven central wavelengths; each config only reports
# the points that fall inside its own fitted range.
SAMPLE = (5800, 6200, 6600, 7000, 7400, 7800)

# The exposure-time floor has to be applied HERE too, not only in run_rm1.py:
# the spec2d files already exist on disk, so the coadd would otherwise sweep the
# short frames back in.  Same yielding rule -- drop them when the configuration
# can spare them, keep them when dropping would leave fewer than two frames and
# trip the single-frame cube bug.
MIN_SCI_EXPTIME = 10.0
MIN_FRAMES_AFTER_CUT = 2


def apply_exptime_floor(spec2ds, tag):
    et = {f: fits.getheader(f).get("EXPTIME", 0.0) for f in spec2ds}
    short = [f for f in spec2ds if et[f] < MIN_SCI_EXPTIME]
    keep = [f for f in spec2ds if f not in short]
    if not short:
        return spec2ds
    if len(keep) >= MIN_FRAMES_AFTER_CUT:
        print(f"  dropping {len(short)} frame(s) under {MIN_SCI_EXPTIME:.0f}s "
              f"({', '.join('%.0fs' % et[f] for f in short)}), "
              f"{len(keep)} remain", flush=True)
        return keep
    print(f"  KEEPING {len(short)} frame(s) under {MIN_SCI_EXPTIME:.0f}s -- "
          f"dropping them would leave {len(keep)}, under the "
          f"{MIN_FRAMES_AFTER_CUT} needed for a valid cube", flush=True)
    return spec2ds


def groups(spec2ds, max_dam, gap_min=None):
    """Split frames into cubes: by STAR first, then by airmass and pointing.

    RM1 is one star throughout, so the star split never fires.  It is kept
    because a driver that silently coadds two different standards is the failure
    RH4 documents, and nothing guarantees the next RM1 night is single-star.

    Within one star: airmass sets the extinction correction, and pypeit applies
    that correction once per cube from a single header airmass (coadd3d.py
    disables the per-frame correction outright).  Measured against PypeIt's
    Mauna Kea curve, an 0.10 airmass spread is ~0.8% in flux, under the ~1%
    scatter this measurement resolves.
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


def star_peak(wl, rad):
    """Spaxel of the STAR in a white-light image, not merely its brightest one.

    `np.argmax` picked the wrong spaxel in two of RM1's nine configs.  2023-12-11
    C peaks on a resampling artifact in row 0, and 2024-03-15 B on a cosmic-ray
    blob 10.6 spaxels (14.4") from the star; both then reported a capture
    fraction measured around the wrong centre -- 67.1% and 40.9% -- and the
    second tripped a "star is not compact" warning about a cube that was fine.
    Neither hurt the science, because pypeit_extract_datacube locates its own
    object, but a guard that cries wolf on good data stops being read.

    Two properties separate a star from both impostors: it is not on the border,
    and it is EXTENDED.  Smoothing by a fraction of the aperture before taking
    the maximum keeps the star, which spans the PSF, and suppresses one- and
    two-spaxel spikes, whose flux is divided over the whole box.
    """
    box = max(2, int(round(rad / 2.0)))
    sm = ndimage.uniform_filter(wl, size=box, mode="nearest")

    # Only positions whose FULL aperture fits inside the cube are candidates.
    # A one-spaxel border is not enough: on 2023-12-11 C the row-0 artifact is
    # 4243 against the star's 1716, so excluding row 0 alone just moves the peak
    # onto the artifact's own smoothed wing at (1,8).  Requiring the aperture to
    # fit is the same condition under which `core` is meaningful at all -- a peak
    # closer to the edge than `rad` has part of its aperture off the cube, so its
    # capture fraction could not be measured even if it were the star.
    m = int(np.ceil(rad))
    ok = np.zeros(wl.shape, bool)
    if 2 * m < min(wl.shape):
        ok[m:-m, m:-m] = True
    else:                           # cube too small for a full aperture inland
        ok[1:-1, 1:-1] = True
    if not ok.any():
        return np.unravel_index(np.argmax(sm), sm.shape)
    return np.unravel_index(np.argmax(np.where(ok, sm, -np.inf)), sm.shape)


def check_cube(cubepath, tag, boxcar=None):
    """Cube geometry guard, carried over from run_rh3_throughput.py unchanged.

    Rejects only a degree-scale CDELT -- the documented single-frame failure that
    inflates extracted counts by ~1e21 and ruins the sensfunc silently.  An
    oversized spatial grid only warns: three RH2 configs came out 553-1307
    spaxels on one axis and were measured to be fine.

    Concentration is measured inside the EXTRACTION APERTURE converted to
    spaxels using THIS cube's own CDELT, not a fixed number of spaxels.  RM1 is
    the run where that matters most: at 3.4" the aperture is 2.5 spaxels on
    Large/2x2 and 10.0 on Small/1x1, so any fixed-spaxel radius would measure
    the star on one night and the field edge on the other.
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

        core, ring, rad, local = np.nan, np.nan, np.nan, np.nan
        wl = np.nansum(np.where(np.isfinite(cd3), cd3, 0), axis=0)
        tot = wl.sum()
        sp = cd[1] * 3600.0 if np.isfinite(cd[1]) else 1.0
        if tot > 0:
            rad = (boxcar / sp) if (boxcar and sp > 0) else 2.5
            j = star_peak(wl, rad)
            yy, xx = np.mgrid[0:wl.shape[0], 0:wl.shape[1]]
            r = np.hypot(yy - j[0], xx - j[1])
            core = wl[r <= rad].sum() / tot
            # Compactness must be judged LOCALLY.  `core` divides by the whole
            # cube, so a bright artifact anywhere inflates the denominator and
            # makes a perfectly compact star look smeared: with the peak-finding
            # fixed, 2023-12-11 C reads 44.1% and 2024-03-15 B 56.0%, both under
            # the 0.60 warning threshold, and both stars are fine.  Against the
            # flux within 3 apertures instead, every RM1 config sits at 93-99%.
            surr = wl[r <= 3 * rad].sum()
            local = (wl[r <= rad].sum() / surr) if surr > 0 else np.nan
            edge = np.zeros(wl.shape, bool)
            edge[0, :] = edge[-1, :] = True
            edge[:, 0] = edge[:, -1] = True
            ring = wl[edge].sum() / tot
        del cd3, wl

        print(f"  cube {shape}  CDELT1/2 = {cd[0]:.3g}/{cd[1]:.3g} deg  "
              f"field {ext[0]:.0f}\" x {ext[1]:.0f}\"  spaxel {sp:.3f}\"  "
              f"flux within {rad:.1f} spaxels ({boxcar}\") of peak "
              f"{100*core:.1f}% of cube, {100*local:.1f}% of surround  "
              f"(outer ring {100*ring:.1f}%)"
              + ("   *** WCS SCALE BROKEN ***" if bad_wcs else ""), flush=True)
        if bad_wcs:
            print("  refusing to extract: a degree-scale CDELT makes counts "
                  "meaningless", flush=True)
            return False
        if np.isfinite(local) and local < 0.60:
            print(f"  WARNING star is not compact ({100*local:.1f}% of the "
                  f"surrounding flux inside the aperture) -- check the "
                  f"whitelight for a smeared or doubled stack", flush=True)
        if np.isfinite(local) and local > 1.05:
            print(f"  NOTE local concentration {100*local:.1f}% exceeds 100%: "
                  f"the annulus around the star carries net NEGATIVE flux, the "
                  f"bowl that accompanies a bright resampling artifact",
                  flush=True)
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
    ap.add_argument("--polyorder", type=int, default=5,
                    help="zeropoint polynomial order.  NOT inherited blindly: "
                         "RL tuned 15 over ~3400 A (227 A per degree of "
                         "freedom) and RH3 showed 15 over ~750 A (47 A/DOF) "
                         "invents shape that independent nights of the same "
                         "star do not share.  SETTLED for RM1 by "
                         "order_study_rm1.py on the 6630 and 7510 night-pairs: "
                         "5 is the only order where BOTH cenwaves hold shape "
                         "agreement under 1.1% with a well-located peak (0.74%/"
                         "58 A at 6630, 1.07%/24 A at 7510), and its ~215 A/DOF "
                         "sits on RL's working 227.  Order 3 wins at 6630 alone "
                         "but scatters 7510's peak by 249 A.  Do not read "
                         "interior scatter to pick this: it falls monotonically "
                         "to order 15 (0.042->0.027) while night-to-night shape "
                         "RMS blows out to 4.9-6.9%, which is what overfitting "
                         "looks like from the inside.")
    ap.add_argument("--boxcar", type=float, default=3.4,
                    help="extraction radius in ARCSEC.  Fixed across cubes on "
                         "purpose: the default aperture is 4*sigma with sigma "
                         "fit in spaxels, and a spaxel is 1.358\" on Large/2x2 "
                         "against 0.339\" on Small/1x1, both of which RM1 uses")
    ap.add_argument("--bridge", type=int, default=5,
                    help="join above-threshold runs separated by at most this "
                         "many pixels before taking the largest.  On RH3 a "
                         "SINGLE pixel reading -184.7 counts split a spectrum "
                         "and the blue half won on length, discarding 353 A of "
                         "good data with nothing in the log to say so")
    ap.add_argument("--suffix", default=None,
                    help="extra tag on the sens filename, for order studies")
    ap.add_argument("--only", default=None)
    ap.add_argument("--skip-coadd", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for sdir in sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run",
                                              "keck_kcrm_*"))):
        if args.only and args.only not in sdir:
            continue
        s2d = sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        if not s2d:
            continue
        s2d = apply_exptime_floor(s2d, sdir)
        night = os.path.relpath(sdir, BASE).split(os.sep)[0]
        letter = os.path.basename(sdir).replace("keck_kcrm_", "")
        log = os.path.join(sdir, "throughput.log")

        # Central wavelength belongs in the tag: 2023-12-10 and 2023-12-11 each
        # carry TWO configurations, so a tag keyed on night alone would collide.
        pf = sorted(glob.glob(os.path.join(sdir, "*.pypeit")))
        m = re.search(r"cenwave:\s*([0-9.]+)", open(pf[0]).read()) if pf else None
        cw = int(round(float(m.group(1)) / 10.0) * 10) if m else 0

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
            print(f"\n=== {tag}  cenwave {cw}: {len(grp)} frames, airmass "
                  f"{am.min():.2f}-{am.max():.2f} (spread {am.max()-am.min():.2f}), "
                  f"span {(mjd.max()-mjd.min())*1440:.0f} min", flush=True)

            cube = f"{tag}.fits"
            c3d = os.path.join(sdir, f"{tag}.coadd3d")
            with open(c3d, "w") as fh:
                fh.write(f"# {star}, {night} setup {letter}, cenwave {cw}; "
                         f"airmass {am.min():.2f}-{am.max():.2f}\n"
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

            sfx = f"p{args.polyorder}cov" + (args.suffix or "")
            sfile = os.path.join(sdir, f"{tag}_IR_{sfx}.sens")
            with open(sfile, "w") as fh:
                fh.write("# IR sensfunc; trim = counts criterion AND full\n"
                         "# slice-coverage interval (measure_trim.py --spec2d).\n"
                         "[sensfunc]\n    algorithm = IR\n    extr = BOX\n"
                         f"    trim_std_pixs = {blue}, {red}\n"
                         f"    polyorder = {args.polyorder}\n")
            sens = f"{tag}_sens_IR_{sfx}.fits"
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
    by magnitudes.  Interior/edge windows are PROPORTIONAL, not a fixed 300 A:
    RH2 showed a fixed window leaves ~90 A on a high-dispersion config and
    reports a meaningless number.  THROUGHPUT is ImageHDU 5 on the WAVE grid,
    NOT a column of the SENS table and NOT on the SENS_WAVE grid.
    """
    with fits.open(path) as hdu:
        d = hdu["SENS"].data
        w = np.asarray(d["SENS_WAVE"]).ravel()
        zp = np.asarray(d["SENS_ZEROPOINT"]).ravel()
        zf = np.asarray(d["SENS_ZEROPOINT_FIT"]).ravel()
        gpm = np.asarray(d["SENS_ZEROPOINT_FIT_GPM"]).ravel().astype(bool)
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
    print(f"  fit range {lo:.0f}-{hi:.0f} A ({hi-lo:.0f} A)   interior scatter "
          f"{sc:.4f} mag   worst edge error {ee:.3f} mag", flush=True)
    if thr is not None and twave is not None and thr.size == twave.size:
        m = (twave >= lo) & (twave <= hi) & np.isfinite(thr) & (thr > 0)
        if m.any():
            tw, tv = twave[m], thr[m]
            q = 0.25 * (tw.max() - tw.min())
            inner = (tw > tw.min() + q) & (tw < tw.max() - q)
            pk = tw[inner][np.argmax(tv[inner])] if inner.any() else np.nan
            print(f"  throughput interior mean {100*np.mean(tv[inner]):.1f}%   "
                  f"peak {100*np.max(tv[inner]):.1f}% @ {pk:.0f} A   "
                  + "  ".join(f"{p}A {100*tv[np.argmin(np.abs(tw-p))]:.1f}%"
                              for p in SAMPLE if lo < p < hi),
                  flush=True)


if __name__ == "__main__":
    main()
