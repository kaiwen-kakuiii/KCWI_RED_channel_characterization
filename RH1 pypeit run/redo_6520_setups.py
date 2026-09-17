#!/usr/bin/env python
"""Re-reduce the two cenwave-6520 setups whose wavelength calibration has
shifted slices, using a one-exposure-wide template.

2023-10-17 C and 2023-10-20 B each ended up with two slices whose solution sits
201-337 A away from their neighbours.  Both were calibrated against
keck_kcrm_RH1.fits, which spans 1103 A -- 1.8 exposures.  full_template cuts a
data-width window out of the template at a position chosen by one
cross-correlation, so a template wider than the data leaves ~1725 px of freedom
in where to cut and dilutes the correlation to cc~0.24.  On those two setups the
correlation preferred a window ~1090 px away for two slices each; the wrong
choices even scored slightly higher cc than the right ones, so nothing flagged
them.

keck_kcrm_RH1_6520.fits is 2064 px spanning 6224-6819 against data at 6216-6820:
template width ~= data width, which leaves essentially no freedom in where to
cut.  That is why every template PypeIt ships is exactly one exposure wide.

Calibrations are deleted so they cannot be reused (PypeIt keys them by frame,
not by parameters).  Science is re-run with -o over the old spec2d.  The
wavelength calibration is checked for agreement between slices BEFORE spending
hours on science -- if a slice is still shifted, the setup stops there.

Usage:
    python redo_6520_setups.py [--dry-run]
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
TEMPLATE = os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH1_6520.fits")
CHECK = os.path.join(ROOT, "pypeit_test", "check_slice_wavelength_agreement.py")

SETUPS = [("2023-10-17", "keck_kcrm_C"), ("2023-10-20", "keck_kcrm_B")]


def sh(cmd, log, cwd):
    with open(log, "a") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(TEMPLATE):
        sys.exit(f"missing {TEMPLATE}")

    for night, setup in SETUPS:
        sdir = os.path.join(BASE, night, "pypeit_run", setup)
        pfile = os.path.join(sdir, f"{setup}.pypeit")
        log = os.path.join(sdir, "redo.log")
        print(f"\n=== {night} {setup}", flush=True)

        txt = open(pfile).read()
        new = re.sub(r"(reid_arxiv = )\S+", r"\1" + TEMPLATE, txt)
        if new == txt and TEMPLATE not in txt:
            print("  could not find reid_arxiv to repoint", flush=True)
            continue
        if args.dry_run:
            print(f"  would point at {os.path.basename(TEMPLATE)} and "
                  f"rebuild {len(glob.glob(os.path.join(sdir, 'Science', 'spec2d_*.fits')))} spec2d")
            continue
        open(pfile, "w").write(new)
        shutil.rmtree(os.path.join(sdir, "Calibrations"), ignore_errors=True)
        print(f"  repointed at {os.path.basename(TEMPLATE)}, cleared Calibrations",
              flush=True)

        if sh([os.path.join(ENV, "run_pypeit"), f"{setup}.pypeit", "-c"], log, sdir):
            print(f"  CALIB FAILED -- see {log}", flush=True)
            continue

        wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
        out = subprocess.run([os.path.join(ENV, "python"), CHECK, "--tol", "50",
                              BASE], capture_output=True, text=True)
        open(log, "a").write(out.stdout + out.stderr)
        line = [l for l in out.stdout.splitlines() if f"{night}/{setup}" in l]
        print(f"  {line[0].strip() if line else 'agreement check produced no line'}",
              flush=True)
        if line and "SHIFTED" in line[0]:
            print("  still shifted -- stopping before science", flush=True)
            continue

        print("  science ...", flush=True)
        rc = sh([os.path.join(ENV, "run_pypeit"), f"{setup}.pypeit", "-o"], log, sdir)
        n = len(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        print(f"  {'DONE' if rc == 0 else f'FAILED rc={rc}'} ({n} spec2d)", flush=True)

    print("\nREDO COMPLETE", flush=True)


if __name__ == "__main__":
    main()
