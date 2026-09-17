#!/usr/bin/env python
"""Decide whether an RH3 template is fit to reduce with, before spending nights on it.

Runs calibrations for one setup with `method = full_template` against the
cenwave-matched template built by `build_rh3_template.py`, and reports how many
of the 24 slices got a usable wavelength solution.  That count is the verdict:
the unsolved slices are what crash the flat field, so a template that leaves any
of them is not usable.

**The template and the line list are one unit, exactly as on RH4.**  An earlier
version of this gate removed `lamps` on the reasoning that RH3's default
FeI/ArI/ArII catalogue is adequate in band (23 lines against 32 detected).  That
reasoning applies to a FeAr arc, and RH3 no longer has one: `seed_rh3.py` retypes
the ThAr frames as the arc and comments the FeAr frames out, because FeAr's 23
lines are too few for holy-grail to identify a solution from scratch.  With ThAr
frames and a ThAr-derived template, the catalogue must be ThAr too.

Measured, on 2023-11-12, which is what exposed it:

    lamps = FeI, ArI, ArII (default)   9-12 lines fitted per slit
    lamps = ThArRH3                    58-78 lines fitted per slit

The default run still reported PASS 24/24 at rms 0.029 px -- low RMS through
about one line per 80 A, which is the "few points, small residual" pattern
RH1_PROCEDURE.md warns about, not a well-constrained solution.  Only the argon
lines of the ThAr spectrum appear in PypeIt's list, so eleven of them carried an
order-4 polynomial across 896 A.

So this gate now sets BOTH, and refuses to run without the line list.

The pre-existing `Calibrations/` is deleted first.  PypeIt keys calibrations by
frame, not by parameters, so the holy-grail products would otherwise be reused
and the test would measure nothing.

## What this gate CAN claim, unlike RH4's

RH2's gate was stronger than RH1's because it ran the template against a
*different night* from the one its seeds came from, which demonstrates that the
template generalises rather than reproducing its own arc.  RH4 could not do that:
one night, one configuration, so its gate re-solved the arc it was seeded from
and could only ever be a self-consistency check.

**RH3 at cenwave 8600 has three independent nights** -- 2023-11-12, 2024-10-29
and 2024-11-06 -- so gating on a night the seeds did not come from is available
here, and `--cross-check` reports which case you are in rather than leaving the
reader to work it out.  Cenwave 8400 has only the one night (2024-10-01) and is
therefore self-consistent, exactly as RH4 was; the gate says so out loud.

Usage:
    python gate_rh3_template.py --night 2024-11-06 --cenwave 8600
    python gate_rh3_template.py --night 2024-10-01 --cenwave 8400
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank, health              # noqa: E402

ENV = "/opt/miniconda3/envs/pypeit/bin"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RH3 pypeit run")

TEMPLATE = {8600: os.path.join(HERE, "keck_kcrm_RH3_8600.fits"),
            8400: os.path.join(HERE, "keck_kcrm_RH3_8400.fits")}
LAMPS = os.path.join(HERE, "ThArRH3")          # PypeIt appends _lines.dat

# Which nights sit at which central wavelength, so the gate can say whether it is
# testing generalisation or reproduction.
NIGHTS = {8600: ["2023-11-12", "2024-10-29", "2024-11-06"],
          8400: ["2024-10-01"]}


def seeded_from(tmpl):
    """Nights whose slices actually seeded this template, read from its header.

    Read from the FILE, not inferred from which WaveCalibs currently look good.
    The inferred version was wrong in a way that flattered the result: once a
    night has been gated, its WaveCalib is a full_template solution and passes
    rank() too, so every night looks like a seed source.  That reported RH3's
    2024-10-29 gate as "independent" when the template was built from precisely
    that night's slices -- the one gate of the three that is reproduction rather
    than generalisation.

    Returns None when the header carries no provenance, so the caller can say
    "unknown" instead of asserting independence it cannot support.
    """
    try:
        with fits.open(tmpl) as h:
            for hdu in h:
                if "SEEDNITE" in hdu.header:
                    return [n for n in str(hdu.header["SEEDNITE"]).split(",") if n]
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", required=True)
    ap.add_argument("--setup", default="keck_kcrm_B")
    ap.add_argument("--cenwave", required=True, type=int, choices=sorted(TEMPLATE))
    ap.add_argument("--allow-cenwave-mismatch", action="store_true",
                    help="gate a setup against a template built for a different "
                         "central wavelength.  Off by default because that is "
                         "precisely why PypeIt's shipped RH3 template fails; "
                         "used deliberately for 2024-10-01, whose 2.3 s arcs are "
                         "too faint to seed a template of its own")
    ap.add_argument("--min-frac", type=float, default=1.0,
                    help="fraction of slices that must solve to pass")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=None,
                    help="physically possible span of one RH3 slice, in A.  "
                         "Measured at 835 A (0.4045 A/px x 2064 px), inside "
                         "pick_seed_slits' 400-900 A default, so this is "
                         "normally unnecessary")
    ap.add_argument("--shift-tol", type=float, default=None,
                    help="A a slice's blue end may sit from the median before "
                         "it counts as SHIFTED.  RH3's genuine slice-to-slice "
                         "spread is the largest measured in this project "
                         "(44.6 A on the 2024-10-29 seeds, against 31.8-38.4 A "
                         "for RH2 and RH4), so the 50 A default has only ~5 A "
                         "of headroom and can reject a correct template")
    args = ap.parse_args()

    tmpl = TEMPLATE[args.cenwave]
    sdir = os.path.join(RUNS, args.night, "pypeit_run", args.setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pfile[0]
    txt = open(pfile).read()

    if not os.path.exists(tmpl):
        sys.exit(f"template missing: {tmpl} -- run "
                 f"build_rh3_template.py --cenwave {args.cenwave} first")
    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"line list missing: {LAMPS}_lines.dat -- the arc frames are "
                 f"ThAr, so the catalogue must be too.  Run:\n"
                 f"  python pypeit_test/build_rh4_linelist.py --label RH3 "
                 f"--wmin 7900 --wmax 9100 --amp-min 300 "
                 f"--out pypeit_test/ThArRH3_lines.dat")

    mc = re.search(r"cenwave:\s*([0-9.]+)", txt)
    setup_cw = int(round(float(mc.group(1)) / 50.0) * 50) if mc else None
    if setup_cw is not None and setup_cw != args.cenwave:
        if not args.allow_cenwave_mismatch:
            sys.exit(f"{args.night}/{args.setup} is at cenwave {setup_cw}, not "
                     f"{args.cenwave} -- gating it against the wrong template "
                     f"would repeat the exact failure this template exists to "
                     f"fix.  Pass --allow-cenwave-mismatch to do it deliberately.")
        # Deliberate, so say exactly how much of the detector has no template
        # under it.  RH3 runs 835 A per slice; the shipped keck_kcrm_RH3.fits
        # failed at ~50% coverage, which is the yardstick for reading this.
        import numpy as _np
        from astropy.io import fits as _fits
        tw = _np.asarray(_fits.getdata(tmpl)["wave"], float)
        lo, hi = setup_cw - 417.0, setup_cw + 417.0
        cov = (min(hi, tw.max()) - max(lo, tw.min())) / (hi - lo)
        print(f"*** DELIBERATE cenwave mismatch: setup is {setup_cw}, template "
              f"is {args.cenwave}")
        print(f"    setup covers ~{lo:.0f}-{hi:.0f} A, template covers "
              f"{tw.min():.0f}-{tw.max():.0f} A -> {100*cov:.0f}% of the "
              f"detector has template under it")
        print(f"    (the shipped keck_kcrm_RH3.fits failed at ~50%)")

    md = re.search(r"decker:\s*(\S+)", txt)
    print(f"{args.night}/{args.setup}: decker {md.group(1) if md else '?'}, "
          f"cenwave {args.cenwave}")
    print(f"   template {os.path.basename(tmpl)}")
    print(f"   lamps    {os.path.basename(LAMPS)}_lines.dat -- the arc frames "
          f"are ThAr, so the catalogue must be ThAr")

    # Rewrite the whole wavelengths block rather than patching pieces of it, so a
    # `lamps` or `method` line left over from the seed stage cannot survive into
    # the gated settings.
    block = ("[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
             "    [[wavelengths]]\n        method = full_template\n"
             f"        reid_arxiv = {tmpl}\n"
             f"        lamps = {LAMPS}\n")
    body = re.sub(r"\[calibrations\]\n    \[\[wavelengths\]\]\n"
                  r"(        \S.*\n|        #.*\n)+", "", txt)
    body = body.replace("[rdx]\n    spectrograph = keck_kcrm\n", block, 1)
    open(pfile, "w").write(body)
    print(f"patched {os.path.basename(pfile)}")

    calib = os.path.join(sdir, "Calibrations")
    if os.path.isdir(calib):
        shutil.rmtree(calib)
        print("removed the holy-grail Calibrations/ so it cannot be reused")

    log = os.path.join(sdir, "gate.log")
    print(f"running calibrations ... (log: {log})", flush=True)
    with open(log, "w") as f:
        rc = subprocess.run([os.path.join(ENV, "run_pypeit"),
                             os.path.basename(pfile), "-c"],
                            stdout=f, stderr=subprocess.STDOUT, cwd=sdir).returncode

    wc = glob.glob(os.path.join(calib, "WaveCalib_*.fits"))
    if not wc:
        sys.exit(f"FAIL: run_pypeit rc={rc} and no WaveCalib was written -- see {log}")

    hkw = {} if args.shift_tol is None else {'shift_tol': args.shift_tol}
    n, solved, flagged, rms, shifted = health(wc[0], span=args.span, **hkw)
    ok = [s for s in solved if s not in flagged and s not in shifted]
    good = rank(wc[0], quiet=True, span=args.span)

    print(f"\nrun_pypeit -c returned {rc}")
    print(f"slices solved, unflagged and unshifted: {len(ok)} / {n}")
    if flagged:
        print(f"   {len(flagged)} flagged by PypeIt: {flagged}")
    if shifted:
        print(f"   {len(shifted)} SHIFTED vs their neighbours: {shifted}"
              f"   (tolerance {args.shift_tol or 50.0:.0f} A)")
    if rms:
        print(f"   rms  median {np.median(rms):.3f}  max {max(rms):.3f} binned px")
    if good:
        print(f"   {len(good)} meet the stricter seed criteria, spanning "
              f"{min(g['wmin'] for g in good):.1f}-{max(g['wmax'] for g in good):.1f} A")

    passed = rc == 0 and len(ok) >= args.min_frac * n
    print(f"\n{'PASS' if passed else 'FAIL'}: template solved {len(ok)}/{n} slices"
          + ("" if passed else f", and run_pypeit returned {rc}"))

    # State plainly which kind of evidence this is, rather than leaving a reader
    # to infer it from the night list.
    seeds = seeded_from(tmpl)
    if passed:
        if seeds is None:
            print(f"NOTE: {os.path.basename(tmpl)} carries no SEEDNITE header, so "
                  f"whether this gate is independent CANNOT be determined -- "
                  f"rebuild it with build_rh3_template.py to stamp provenance.")
        elif args.night in seeds:
            print(f"NOTE: self-consistency only -- {args.night} is one of the "
                  f"nights this template was seeded from ({', '.join(seeds)}), so "
                  f"this re-solves the arc it came from, as RH4's gate did.")
        else:
            print(f"NOTE: INDEPENDENT -- this template was seeded from "
                  f"{', '.join(seeds)}, not {args.night}, so passing here "
                  f"demonstrates generalisation rather than reproduction.")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
