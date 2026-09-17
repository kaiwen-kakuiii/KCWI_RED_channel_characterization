#!/usr/bin/env python
"""Drive the RH1 standard-star reduction, night by night, with the RL procedure.

The RH1 standards are g191b2b (7 nights) and feige34 (2 nights), not the
feige110 the RL run used; nothing here depends on which, since the star is
located per setup rather than assumed.

Per setup:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RH1/<night> -c all
  2. for each setup that contains science frames:
       a. run_pypeit -c                      (calibrations only -> Slits)
       b. find_object_regions.py            (star position -> user_regions)
       c. insert [reduce][[skysub]] user_regions into the .pypeit
       d. run_pypeit                        (science, reusing calibs)

Sequential on purpose: one run_pypeit saturates the machine, and parallel
sensfuncs already bit us once (telluric cache race, see THROUGHPUT_PROCEDURE).

Idempotent-ish: a setup whose Science/ already has spec2d files is skipped,
so a crashed driver can be relaunched.
"""
import os
import re
import glob
import shutil
import subprocess
import sys

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
RAW = os.path.join(ROOT, "fits", "by_night", "RH1")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
TEMPLATE = os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH1.fits")

def sh(cmd, log, cwd=None):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=cwd).returncode

def science_rows(pfile):
    rows = []
    for line in open(pfile):
        m = re.match(r"\s*(\S+\.fits)\s*\|\s*([\w,]+)\s*\|", line)
        if m and "science" in m.group(2).split(","):
            rows.append(m.group(1))
    return rows

for night in sorted(os.listdir(RAW)):
    ndir = os.path.join(BASE, night)
    os.makedirs(ndir, exist_ok=True)
    log = os.path.join(ndir, "driver.log")
    rundir = os.path.join(ndir, "pypeit_run")
    if not glob.glob(os.path.join(rundir, "keck_kcrm_*", "*.pypeit")):
        rc = sh([os.path.join(ENV, "pypeit_setup"), "-s", "keck_kcrm",
                 "-r", os.path.join(RAW, night), "-c", "all", "-d", rundir], log)
        if rc:
            print(f"{night}: SETUP FAILED rc={rc}", flush=True)
            continue

    for sdir in sorted(glob.glob(os.path.join(rundir, "keck_kcrm_*"))):
        pfile = glob.glob(os.path.join(sdir, "*.pypeit"))[0]
        sci = science_rows(pfile)
        tag = f"{night}/{os.path.basename(sdir)}"
        if not sci:
            print(f"{tag}: no science frames, skipped", flush=True)
            continue
        if glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")):
            print(f"{tag}: spec2d already present, skipped", flush=True)
            continue

        # PypeIt ships no RH1 template, so without this it silently falls back to
        # holy-grail and solves ~1 of 24 slices -- see build_rh1_template.py.
        # Calibrations are keyed by frame, not by parameters, so a Calibrations/
        # built before this block would be *reused* rather than redone: any
        # directory predating the template is thrown away instead.
        # Test for *any* reid_arxiv, not for this template's filename: the Small
        # slicer setup is pointed at keck_kcrm_RH1_small.fits, and a check for the
        # default name would not match it -- so this block would repatch a setup
        # that is already configured and delete calibrations that took hours.
        txt = open(pfile).read()
        if "reid_arxiv" not in txt:
            txt = txt.replace(
                "[rdx]\n    spectrograph = keck_kcrm",
                "[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
                "    [[wavelengths]]\n        method = full_template\n"
                f"        reid_arxiv = {TEMPLATE}", 1)
            open(pfile, "w").write(txt)
            if os.path.isdir(os.path.join(sdir, "Calibrations")):
                shutil.rmtree(os.path.join(sdir, "Calibrations"))
                print(f"{tag}: cleared pre-template Calibrations", flush=True)

        print(f"{tag}: {len(sci)} science frames -- calibrations ...", flush=True)
        if sh([os.path.join(ENV, "run_pypeit"), os.path.basename(pfile), "-c"],
              log, cwd=sdir):
            print(f"{tag}: CALIB RUN FAILED, night continues with next setup", flush=True)
            continue

        # How many of the 24 slices actually got a wavelength solution.  A run
        # that gets this far has not crashed, but a setup reduced with a third of
        # its slices masked is not the same measurement as one with all of them,
        # and nothing downstream says so out loud.
        wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
        if wc:
            chk = subprocess.run([os.path.join(ENV, "python"), SEEDS, wc[0]],
                                 capture_output=True, text=True)
            open(log, "a").write(chk.stdout + chk.stderr)
            m = re.search(r"(\d+) usable", chk.stdout)
            print(f"{tag}: wavelength solutions on "
                  f"{m.group(1) if m else '?'} slices", flush=True)

        # star position -> sky regions, from the first science frame
        out = subprocess.run(
            [os.path.join(ENV, "python"), FINDER,
             os.path.join(RAW, night, sci[0]), os.path.join(sdir, "Calibrations")],
            capture_output=True, text=True)
        open(log, "a").write(out.stdout + out.stderr)
        m = re.search(r"user_regions = (\S+)", out.stdout)
        if not m:
            print(f"{tag}: FIND_OBJECT_REGIONS FAILED, setup skipped", flush=True)
            continue
        regions = m.group(1)
        print(f"{tag}: user_regions = {regions}", flush=True)

        txt = open(pfile).read()
        if "[[skysub]]" not in txt:
            txt = txt.replace(
                f"        reid_arxiv = {TEMPLATE}",
                f"        reid_arxiv = {TEMPLATE}\n\n[reduce]\n"
                f"    [[skysub]]\n        user_regions = {regions}", 1)
            open(pfile, "w").write(txt)

        print(f"{tag}: science run ...", flush=True)
        rc = sh([os.path.join(ENV, "run_pypeit"), os.path.basename(pfile)],
                log, cwd=sdir)
        n2d = len(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        print(f"{tag}: {'DONE' if rc == 0 else f'FAILED rc={rc}'} "
              f"({n2d}/{len(sci)} spec2d)", flush=True)

print("ALL NIGHTS PROCESSED", flush=True)
