#!/usr/bin/env python
"""Decide whether the RH4 template is fit to reduce with, before spending the night on it.

Runs calibrations for one setup with `method = full_template` against
`keck_kcrm_RH4.fits` **and** the ThArRH4 line list, and reports how many of the
24 slices got a usable wavelength solution.  That count is the verdict: the
unsolved slices are what crash the flat field, so a template that leaves any of
them is not usable.

The pre-existing `Calibrations/` is deleted first.  PypeIt keys calibrations by
frame, not by parameters, so the holy-grail products would otherwise be reused
and the test would measure nothing.

**Both settings are checked, not just the template.**  RH4's solution is a
template *and* a line list; PypeIt's KCWI default catalogue has four lines above
9550 A, so a run that quietly kept `lamps = FeI, ArI, ArII` would fail in a way
that looks like a bad template.

## What this gate cannot do

RH2's gate was stronger than RH1's because it ran the template against a
*different night* from the one its seeds came from, which demonstrates that the
template generalises rather than reproducing its own arc.

**RH4 has one night and one configuration, so that test does not exist here.**
This gate re-solves the same arc the template was seeded from.  It still catches
the failure it is meant to catch -- holy-grail solving a handful of slices and
the flat field dying on the rest -- and the improvement it measures is real,
because identification changes while the fit is to catalogue wavelengths either
way (RH1_WAVELENGTH.md).  But it is a self-consistency check, not a
generalisation test, and the reduction should be read with that in mind.  If a
second RH4 night is ever fetched, gate on it.

Usage:
    python gate_rh4_template.py --night 2024-12-28 --setup keck_kcrm_B
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank, health              # noqa: E402

ENV = "/opt/miniconda3/envs/pypeit/bin"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RH4 pypeit run")

TEMPLATE = os.path.join(HERE, "keck_kcrm_RH4.fits")
LAMPS = os.path.join(HERE, "ThArRH4")          # PypeIt appends _lines.dat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="2024-12-28")
    ap.add_argument("--setup", default="keck_kcrm_B")
    ap.add_argument("--min-frac", type=float, default=1.0,
                    help="fraction of slices that must solve to pass")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=None,
                    help="physically possible span of one RH4 slice, in A")
    args = ap.parse_args()

    sdir = os.path.join(RUNS, args.night, "pypeit_run", args.setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pfile[0]
    txt = open(pfile).read()

    if not os.path.exists(TEMPLATE):
        sys.exit(f"template missing: {TEMPLATE} -- run build_rh4_template.py first")
    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"line list missing: {LAMPS}_lines.dat -- run "
                 f"build_rh4_linelist.py first")

    m = re.search(r"decker:\s*(\S+)", txt)
    print(f"{args.night}/{args.setup}: decker {m.group(1) if m else '?'}")
    print(f"   template {os.path.basename(TEMPLATE)}")
    print(f"   lamps    {os.path.basename(LAMPS)}_lines.dat")

    # Rewrite the whole wavelengths block rather than patching pieces of it: the
    # template and the line list have to agree, and a file left half-repointed
    # would be gated on settings nobody chose.
    block = ("[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
             "    [[wavelengths]]\n        method = full_template\n"
             f"        reid_arxiv = {TEMPLATE}\n"
             f"        lamps = {LAMPS}\n")
    body = re.sub(r"\[calibrations\]\n    \[\[wavelengths\]\]\n(        \S.*\n)+",
                  "", txt)
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

    n, solved, flagged, rms, shifted = health(wc[0], span=args.span)
    ok = [s for s in solved if s not in flagged and s not in shifted]
    good = rank(wc[0], quiet=True, span=args.span)

    print(f"\nrun_pypeit -c returned {rc}")
    print(f"slices solved, unflagged and unshifted: {len(ok)} / {n}")
    if flagged:
        print(f"   {len(flagged)} flagged by PypeIt: {flagged}")
    if shifted:
        print(f"   {len(shifted)} SHIFTED vs their neighbours: {shifted}")
    if rms:
        print(f"   rms  median {np.median(rms):.3f}  max {max(rms):.3f} binned px")
    if good:
        print(f"   {len(good)} meet the stricter seed criteria, spanning "
              f"{min(g['wmin'] for g in good):.1f}-{max(g['wmax'] for g in good):.1f} A")

    passed = rc == 0 and len(ok) >= args.min_frac * n
    print(f"\n{'PASS' if passed else 'FAIL'}: template solved {len(ok)}/{n} slices"
          + ("" if passed else f", and run_pypeit returned {rc}"))
    if passed:
        print("NOTE: this is a self-consistency check -- RH4 has only the one "
              "night, so the template was gated on the arc it was seeded from.")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
