#!/usr/bin/env python
"""Take the RH1 spec2d files through to throughput, per THROUGHPUT_PROCEDURE.md.

    spec2d  --pypeit_coadd_datacube-->  datacube + whitelight
            --pypeit_extract_datacube-> spec1d (BOX aperture; --boxcar ARCSEC)
            --measure_trim.py --spec2d-> trim_std_pixs
            --pypeit_sensfunc---------> zeropoint + throughput

One cube per VISIT, not per setup, split on AIRMASS -- see visits() for why that
is the criterion and how the threshold is derived.  Five of the ten RH1 setups
were observed in two blocks four to five hours apart, spanning airmass 1.21-1.63,
which is 2-3% in flux at the blue end of RH1 and would dominate the ~1% scatter
this measurement resolves.  2023-10-15 must be split regardless of airmass: its
two blocks sit at different pointings (ra_off 3.3 vs 0.0), which `combine = True`
would smear rather than stack.

Sensfuncs run one at a time on purpose: running three wide, one job found a
half-written entry in the telluric cache and died (see THROUGHPUT_PROCEDURE.md).

Usage:
    python run_rh1_throughput.py [--max-dam 0.10] [--gap MIN] [--polyorder 15]
                                 [--boxcar ARCSEC] [--skip-coadd] [--dry-run]

To re-extract every cube through one fixed aperture, cubes already built:

    python run_rh1_throughput.py --boxcar 3.4 --skip-coadd
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


def visits(spec2ds, max_dam, gap_min=None):
    """Group frames so that each cube has a well-defined airmass.

    Airmass is the criterion, not elapsed time: it is what sets the extinction
    correction, and pypeit applies that correction once per cube from a single
    header airmass rather than per frame (coadd3d.py disables the per-frame
    correction outright, noting the standard star exposures "are assumed to have
    similar airmasses").  Time is only a proxy -- frames minutes apart near the
    horizon can differ more in airmass than frames an hour apart near transit.

    Measured against PypeIt's own Mauna Kea extinction curve at 6200 A, the
    bluest RH1 coverage and so the worst case:

        airmass spread   0.05    0.10    0.15    0.20    0.36
        flux error      0.39%   0.78%   1.17%   1.57%   2.84%

    The default 0.10 keeps the induced error near 0.8%, under the ~1% scatter
    this measurement resolves.

    Pointing is a hard split regardless of airmass: frames at different ra_off /
    dec_off would be smeared rather than stacked by `combine = True`.
    """
    rows = []
    for f in spec2ds:
        h = fits.getheader(f)
        rows.append((h.get("MJD", 0.0), h.get("AIRMASS", 0.0),
                     (round(h.get("RAOFF", 0.0) or 0.0, 2),
                      round(h.get("DECOFF", 0.0) or 0.0, 2)), f))
    rows.sort()

    groups, cur = [], [rows[0]]
    for prev, this in zip(rows, rows[1:]):
        ams = [r[1] for r in cur] + [this[1]]
        split = (this[2] != cur[0][2]                       # pointing changed
                 or max(ams) - min(ams) > max_dam)          # airmass drifted
        if gap_min is not None and (this[0] - prev[0]) * 1440.0 > gap_min:
            split = True
        if split:
            groups.append(cur)
            cur = []
        cur.append(this)
    groups.append(cur)
    return groups


def sh(cmd, log, cwd):
    with open(log, "a") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-dam", type=float, default=0.10,
                    help="largest airmass spread allowed inside one cube; 0.10 "
                         "holds the extinction error near 0.8%% at 6200 A")
    ap.add_argument("--gap", type=float, default=None,
                    help="optional extra split on a time gap, in minutes")
    ap.add_argument("--polyorder", type=int, default=15)
    ap.add_argument("--only", default=None,
                    help="substring a setup directory must contain, so a rerun "
                         "can redo just the setups that changed")
    ap.add_argument("--boxcar", type=float, default=None,
                    help="extraction radius in ARCSEC, passed to "
                         "pypeit_extract_datacube -b.  Fix this across nights "
                         "when comparing throughput: the default aperture is "
                         "4*sigma with sigma fit in SPAXELS, and the spaxel is "
                         "1.358 arcsec on Large/2x2 but 0.339 arcsec on "
                         "Small/1x1, so the default measures a different solid "
                         "angle per night.  The fit also floors sigma at 0.5 "
                         "spaxel, a floor that binds only on the coarse grid.")
    ap.add_argument("--skip-coadd", action="store_true",
                    help="reuse the existing datacube instead of rebuilding it; "
                         "for reruns that only change the extraction")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    setups = sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run", "keck_kcrm_*")))
    for sdir in setups:
        if args.only and args.only not in sdir:
            continue
        s2d = sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        if not s2d:
            continue
        night = os.path.relpath(sdir, BASE).split(os.sep)[0]
        letter = os.path.basename(sdir).replace("keck_kcrm_", "")
        star = re.search(r"spec2d_\S+?-(\S+?)_KCRM", os.path.basename(s2d[0])).group(1)
        log = os.path.join(sdir, "throughput.log")

        groups = visits(s2d, args.max_dam, gap_min=args.gap)
        for gi, grp in enumerate(groups, start=1):
            am = np.array([r[1] for r in grp])
            mjd = np.array([r[0] for r in grp])
            tag = f"{star}_{night}_{letter}" + (f"_v{gi}" if len(groups) > 1 else "")
            print(f"\n=== {tag}: {len(grp)} frames, airmass "
                  f"{am.min():.2f}-{am.max():.2f} (spread {am.max()-am.min():.2f}), "
                  f"span {(mjd.max()-mjd.min())*1440:.0f} min", flush=True)

            cube = f"{tag}.fits"
            c3d = os.path.join(sdir, f"{tag}.coadd3d")
            with open(c3d, "w") as fh:
                fh.write(f"# {star}, {night} setup {letter}, visit {gi} of "
                         f"{len(groups)}; airmass {am.min():.2f}-{am.max():.2f}\n"
                         "[rdx]\n    spectrograph = keck_kcrm\n\n[reduce]\n"
                         "    [[cube]]\n        combine = True\n"
                         f"        output_filename = {cube}\n"
                         "        save_whitelight = True\n\nspec2d read\nfilename\n")
                for r in grp:
                    f = r[-1]
                    fh.write(f"Science/{os.path.basename(f)}\n")
                fh.write("spec2d end\n")
            if args.dry_run:
                continue

            if args.skip_coadd and os.path.exists(os.path.join(sdir, cube)):
                print("  reusing existing cube", flush=True)
            elif sh([os.path.join(ENV, "pypeit_coadd_datacube"),
                     os.path.basename(c3d), "-o"], log, sdir):
                print(f"  COADD FAILED -- see {log}", flush=True)
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
            out = subprocess.run(
                [os.path.join(ENV, "python"), TRIM, "--spec2d",
                 os.path.join(sdir, "Science"), spec1d[-1]],
                capture_output=True, text=True)
            open(log, "a").write(out.stdout + out.stderr)
            m = re.search(r"trim_std_pixs = (\d+), (\d+)", out.stdout)
            if not m:
                print(f"  TRIM FAILED -- see {log}", flush=True)
                continue
            blue, red = m.group(1), m.group(2)

            sfile = os.path.join(sdir, f"{tag}_IR_p{args.polyorder}cov.sens")
            with open(sfile, "w") as fh:
                fh.write("# IR sensfunc; trim = counts criterion AND full\n"
                         "# slice-coverage interval (see measure_trim.py --spec2d).\n"
                         "[sensfunc]\n    algorithm = IR\n    extr = BOX\n"
                         f"    trim_std_pixs = {blue}, {red}\n"
                         f"    polyorder = {args.polyorder}\n")
            sens = f"{tag}_sens_IR_p{args.polyorder}cov.fits"
            rc = sh([os.path.join(ENV, "pypeit_sensfunc"),
                     os.path.basename(spec1d[-1]), "-s", os.path.basename(sfile),
                     "-o", sens], log, sdir)
            print(f"  trim {blue},{red}  -> "
                  f"{'sensfunc OK' if rc == 0 else f'SENSFUNC FAILED rc={rc}'}: {sens}",
                  flush=True)

    print("\nALL SETUPS PROCESSED", flush=True)


if __name__ == "__main__":
    main()
