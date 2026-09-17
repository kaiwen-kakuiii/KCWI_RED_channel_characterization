#!/usr/bin/env python
"""Decide whether the stitched RH1 template is fit to reduce nine nights with.

Runs calibrations for one setup with `method = full_template` pointed at the new
template, and reports how many of the 24 slices got a usable wavelength solution.
Holy-grail managed 1 of 24 on this same setup, and the 23 failures are what
crashed the flat field, so this is the number that matters.

The pre-existing Calibrations/ directory is deleted first: PypeIt keys
calibrations by frame, not by parameters, so the holy-grail products would
otherwise be reused and the test would measure nothing.

Usage:
    python gate_rh1_template.py [--night 2023-09-16] [--setup keck_kcrm_B]
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank                       # noqa: E402

ENV = "/opt/miniconda3/envs/pypeit/bin"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TEMPLATE = os.path.join(HERE, "keck_kcrm_RH1.fits")

PARAMS = ("[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
          "    [[wavelengths]]\n        method = full_template\n"
          f"        reid_arxiv = {TEMPLATE}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="2023-09-16")
    ap.add_argument("--setup", default="keck_kcrm_B")
    ap.add_argument("--min-frac", type=float, default=1.0,
                    help="fraction of slices that must solve to pass")
    args = ap.parse_args()

    sdir = os.path.join(ROOT, "RH1 pypeit run", args.night, "pypeit_run", args.setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pfile[0]
    if not os.path.exists(TEMPLATE):
        sys.exit(f"template missing: {TEMPLATE} -- run build_rh1_template.py first")

    txt = open(pfile).read()
    if "keck_kcrm_RH1.fits" not in txt:
        open(pfile, "w").write(
            txt.replace("[rdx]\n    spectrograph = keck_kcrm", PARAMS, 1))
        print(f"patched {os.path.basename(pfile)} to use the template")
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

    good = rank(wc[0], quiet=True)
    from pypeit.wavecalib import WaveCalib
    n = len(WaveCalib.from_file(wc[0], chk_version=False).spat_ids)
    rms = np.array([g["rms"] for g in good])

    print(f"\nrun_pypeit -c returned {rc}")
    print(f"slices with a usable solution: {len(good)} / {n}")
    if len(good):
        print(f"   rms  median {np.median(rms):.3f}  max {rms.max():.3f} binned px")
        print(f"   span {min(g['wmin'] for g in good):.1f}-"
              f"{max(g['wmax'] for g in good):.1f} A")

    passed = rc == 0 and len(good) >= args.min_frac * n
    print(f"\n{'PASS' if passed else 'FAIL'}: "
          f"holy-grail managed 1/{n} on this setup; template managed {len(good)}/{n}"
          + ("" if passed else f", and run_pypeit returned {rc}"))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
