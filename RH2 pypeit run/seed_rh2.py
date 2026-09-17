#!/usr/bin/env python
"""Bootstrap wavelength seeds for RH2, which PypeIt ships no template for.

Same hole as RH1: keck_kcwi.py config_specific_par assigns reid_arxiv for RL,
RM1, RM2 and RH3 only, so RH2 falls back to holy-grail.  Holy-grail solves a
handful of the 24 slices and the unsolved ones then kill the flat field
(flatfield.py illum_profile_spectral takes np.min of an empty array), so each
run here is EXPECTED to end with a non-zero rc -- WaveCalib is written before
the flat field, which is all we need.

RH2 splits by slicer, and the slicer sets spectral resolution, so a template
built at one slicer may not identify lines at another (RH1's Small-slicer setup
failed against a Large-slicer template at a 2.06x line-width ratio).  Seeds are
therefore collected per slicer:

    Large   2023-11-17 B (6900)  +  2023-07-16 A (7180)   -> spans ~6585-7495 A
    Medium  2024-11-07 B (6750)                           -> one exposure wide

Sequential: an RH1 reduction is running alongside this on the same 16 GB box.
"""
import glob
import os
import subprocess
import sys

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")

# (night, setup, slicer, cenwave) -- chosen to span each slicer's range with as
# few holy-grail runs as possible, since each one is slow and mostly fails.
TARGETS = [
    ("2023-11-17", "keck_kcrm_B", "Large",  6900),
    ("2023-07-16", "keck_kcrm_A", "Large",  7180),
    ("2024-11-07", "keck_kcrm_B", "Medium", 6750),
]

for night, setup, slicer, cwave in TARGETS:
    sdir = os.path.join(BASE, night, "pypeit_run", setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        print(f"{night}/{setup}: no .pypeit, skipped", flush=True)
        continue
    pfile = pfile[0]
    wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
    if wc:
        print(f"{night}/{setup}: WaveCalib already present, not rerun", flush=True)
    else:
        print(f"{night}/{setup} ({slicer}, {cwave} A): holy-grail ...", flush=True)
        log = os.path.join(sdir, "holygrail.log")
        with open(log, "w") as f:
            rc = subprocess.run([os.path.join(ENV, "run_pypeit"),
                                 os.path.basename(pfile), "-c"],
                                stdout=f, stderr=subprocess.STDOUT,
                                cwd=sdir).returncode
        wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
        print(f"{night}/{setup}: rc={rc} (non-zero expected), "
              f"WaveCalib {'written' if wc else 'MISSING'}", flush=True)
    if wc:
        out = subprocess.run([os.path.join(ENV, "python"), SEEDS, wc[0]],
                             capture_output=True, text=True)
        print(out.stdout + out.stderr, flush=True)

print("SEED PHASE DONE", flush=True)
