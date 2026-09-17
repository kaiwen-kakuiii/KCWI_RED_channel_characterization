#!/usr/bin/env python
"""Decide whether an RH2 template is fit to reduce with, before spending nights on it.

Runs calibrations for one setup with `method = full_template` pointed at the
template for that setup's slicer, and reports how many of the 24 slices got a
usable wavelength solution.  Holy-grail manages a handful on these setups, and
the unsolved slices are what crash the flat field, so this is the number that
decides.

The pre-existing Calibrations/ is deleted first: PypeIt keys calibrations by
frame, not by parameters, so the holy-grail products would otherwise be reused
and the test would measure nothing.

Usage:
    python gate_rh2_template.py --night 2023-11-17 --setup keck_kcrm_B
    python gate_rh2_template.py --night 2024-11-07 --setup keck_kcrm_B
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
RUNS = os.path.join(ROOT, "RH2 pypeit run")

TEMPLATE = {"Large":  os.path.join(HERE, "keck_kcrm_RH2.fits"),
            "Medium": os.path.join(HERE, "keck_kcrm_RH2_medium.fits")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", required=True)
    ap.add_argument("--setup", default="keck_kcrm_B")
    ap.add_argument("--min-frac", type=float, default=1.0,
                    help="fraction of slices that must solve to pass")
    args = ap.parse_args()

    sdir = os.path.join(RUNS, args.night, "pypeit_run", args.setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pfile[0]
    txt = open(pfile).read()

    # The slicer decides which template applies -- it sets spectral resolution,
    # and a template built at another slicer cannot name lines it does not
    # resolve.  Read it from the setup block rather than assuming.
    m = re.search(r"decker:\s*(\S+)", txt)
    if not m or m.group(1) not in TEMPLATE:
        sys.exit(f"decker {m.group(1) if m else '?'} has no RH2 template defined")
    decker = m.group(1)
    tmpl = TEMPLATE[decker]
    if not os.path.exists(tmpl):
        sys.exit(f"template missing: {tmpl} -- run build_rh2_template.py "
                 f"--decker {decker} first")
    print(f"{args.night}/{args.setup}: decker {decker} -> {os.path.basename(tmpl)}")

    params = ("[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
              "    [[wavelengths]]\n        method = full_template\n"
              f"        reid_arxiv = {tmpl}")
    if "reid_arxiv" not in txt:
        open(pfile, "w").write(
            txt.replace("[rdx]\n    spectrograph = keck_kcrm", params, 1))
        print(f"patched {os.path.basename(pfile)} to use the template")
    elif tmpl not in txt:
        # Pointed at a different template than the one being gated -- repoint it,
        # otherwise the gate reports on a file it did not test.
        open(pfile, "w").write(re.sub(r"reid_arxiv = \S+", f"reid_arxiv = {tmpl}", txt, 1))
        print(f"repointed {os.path.basename(pfile)} at {os.path.basename(tmpl)}")

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

    # health() also returns the slices that sit far from their neighbours in
    # wavelength; a slice can keep a normal width, dispersion and RMS while being
    # matched to the wrong window, and only its neighbours expose that.
    n, solved, flagged, rms, shifted = health(wc[0])
    ok = [s for s in solved if s not in flagged and s not in shifted]
    good = rank(wc[0], quiet=True)

    print(f"\nrun_pypeit -c returned {rc}")
    print(f"slices solved, unflagged and unshifted: {len(ok)} / {n}")
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
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
