#!/usr/bin/env python
"""Drive the RH2 standard-star reduction, night by night, with the RL procedure.

The RH2 standards are feige110 (5 nights) and g191b2b (2023-09-16, added later
at Kaiwen's request to add central wavelengths 6950 and 7100 at the Large
slicer).  Nothing here depends on which star: the object is located per setup by
find_object_regions.py, and pypeit_sensfunc matches the standard on coordinates.

Per setup:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RH2/<night> -c all
  2. for each setup that is COMPLETE and contains science frames:
       a. run_pypeit -c                     (calibrations only -> Slits)
       b. find_object_regions.py            (star position -> user_regions)
       c. insert [reduce][[skysub]] user_regions into the .pypeit
       d. run_pypeit                        (science, reusing calibs)

Three things differ from run_rh1.py, all of them properties of this data set:

**Two templates, chosen by slicer.**  PypeIt ships none for RH2 (the RL/RM1/RM2/
RH3 list in keck_kcwi.py config_specific_par has the same hole that stopped RH1),
and RH2's nights split across two slicers.  The slicer sets spectral resolution,
so a template built at one holds as single features the blends the other
resolves; RH1's Small-slicer setup failed exactly that way against a
Large-slicer template.  Large nights get keck_kcrm_RH2.fits, Medium nights
keck_kcrm_RH2_medium.fits, read from each setup's own decker.

**Incomplete setups are skipped by inspection, not by crashing.**  pypeit_setup
groups on (dispname, decker, binning, cenwave), and a stray frame whose header
disagrees with the rest of its night lands in a setup of its own with no
calibrations -- 2023-11-17 D is one science frame, decker "unknown", seven
biases and nothing else.  Running it wastes an hour to reach a confusing error,
so a setup without both arc and trace frames is reported and passed over.

**Hand-edited .pypeit files are preserved.**  2024-11-06's dome flats carry
FLSPECTR='off' although they are plainly illuminated (24267 ADU against 24038-
24575 on the known-good night), so PypeIt typed them None and left that night
with no trace frames at all.  They were retyped illumflat,trace by hand.  Setup
is therefore only run when no .pypeit exists yet -- regenerating would silently
revert that edit and take the night's only complete config with it.

Sequential on purpose: one run_pypeit saturates the machine, and parallel
sensfuncs already bit us once (telluric cache race, see THROUGHPUT_PROCEDURE).
This driver is meant to run alongside the RH1 one, which is two run_pypeit
processes on a 16 GB box -- not more.

Idempotent-ish: a setup whose Science/ already has spec2d files is skipped, so a
crashed driver can be relaunched.
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
RAW = os.path.join(ROOT, "fits", "by_night", "RH2")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
# Setups the driver refused to reduce, appended as it goes, so a run left
# unattended still leaves a list of what needs a decision.
NEEDS_DECISION = os.path.join(BASE, "NEEDS_DECISION.txt")
TEMPLATE = {"Large":  os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH2.fits"),
            "Medium": os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH2_medium.fits")}


def sh(cmd, log, cwd=None):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=cwd).returncode


# Anything shorter than this is an acquisition or focus test, not a standard
# exposure.  2023-11-07 carried a 5 s frame that sorted *first* among its science
# rows, so it would have been the frame handed to the object finder.
MIN_SCI_EXPTIME = 10.0


def frames(pfile):
    """({frametype: [filename, ...]}, {filename: exptime}) for uncommented rows.

    Columns are located by name from the data block's own header line rather than
    by position: a .pypeit written by another PypeIt version can order them
    differently, and reading exptime out of the wrong column silently drops good
    frames.
    """
    cols, types, exp = None, {}, {}
    for line in open(pfile):
        if line.lstrip().startswith("#"):
            continue
        parts = [q.strip() for q in line.split("|")]
        if cols is None:
            if parts and parts[0] == "filename":
                cols = parts
            continue
        if not parts or not parts[0].endswith(".fits"):
            continue
        row = dict(zip(cols, parts))
        fn = row["filename"]
        for t in row.get("frametype", "").split(","):
            types.setdefault(t, []).append(fn)
        try:
            exp[fn] = float(row.get("exptime", "nan"))
        except ValueError:
            exp[fn] = float("nan")
    return types, exp


for night in sorted(os.listdir(RAW)):
    if not os.path.isdir(os.path.join(RAW, night)):
        continue
    ndir = os.path.join(BASE, night)
    os.makedirs(ndir, exist_ok=True)
    log = os.path.join(ndir, "driver.log")
    rundir = os.path.join(ndir, "pypeit_run")
    # Only when nothing exists: see the note on hand-edited files above.
    if not glob.glob(os.path.join(rundir, "keck_kcrm_*", "*.pypeit")):
        rc = sh([os.path.join(ENV, "pypeit_setup"), "-s", "keck_kcrm",
                 "-r", os.path.join(RAW, night), "-c", "all", "-d", rundir], log)
        if rc:
            print(f"{night}: SETUP FAILED rc={rc}", flush=True)
            continue

    for sdir in sorted(glob.glob(os.path.join(rundir, "keck_kcrm_*"))):
        pfile = glob.glob(os.path.join(sdir, "*.pypeit"))[0]
        tag = f"{night}/{os.path.basename(sdir)}"
        ftypes, exptime = frames(pfile)
        sci = ftypes.get("science", [])
        short = [f for f in sci if not (exptime.get(f, 0.0) >= MIN_SCI_EXPTIME)]
        if short:
            print(f"{tag}: ignoring {len(short)} science frame(s) under "
                  f"{MIN_SCI_EXPTIME:.0f}s -- test frames, not standards", flush=True)
            sci = [f for f in sci if f not in short]
        if not sci:
            print(f"{tag}: no science frames, skipped", flush=True)
            continue
        missing = [t for t in ("arc", "trace", "pixelflat") if not ftypes.get(t)]
        if missing:
            print(f"{tag}: {len(sci)} science but no {'/'.join(missing)} "
                  f"frames -- incomplete config, skipped", flush=True)
            continue
        if glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")):
            print(f"{tag}: spec2d already present, skipped", flush=True)
            continue

        txt = open(pfile).read()
        m = re.search(r"decker:\s*(\S+)", txt)
        decker = m.group(1) if m else None
        if decker not in TEMPLATE:
            print(f"{tag}: decker {decker} has no RH2 template, skipped", flush=True)
            continue
        tmpl = TEMPLATE[decker]
        if not os.path.exists(tmpl):
            print(f"{tag}: {os.path.basename(tmpl)} missing -- run "
                  f"build_rh2_template.py --decker {decker}, skipped", flush=True)
            continue

        # Without this PypeIt falls back to holy-grail silently and solves a
        # handful of the 24 slices, after which the flat field dies on the rest.
        # Test for *any* reid_arxiv rather than this file's name, so a setup
        # already pointed at the other slicer's template is not repatched and its
        # calibrations thrown away.  Calibrations are keyed by frame, not by
        # parameters, so a directory predating the template would be reused.
        if "reid_arxiv" not in txt:
            txt = txt.replace(
                "[rdx]\n    spectrograph = keck_kcrm",
                "[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
                "    [[wavelengths]]\n        method = full_template\n"
                f"        reid_arxiv = {tmpl}", 1)
            open(pfile, "w").write(txt)
            if os.path.isdir(os.path.join(sdir, "Calibrations")):
                shutil.rmtree(os.path.join(sdir, "Calibrations"))
                print(f"{tag}: cleared pre-template Calibrations", flush=True)

        print(f"{tag}: {decker}, {len(sci)} science frames -- calibrations ...",
              flush=True)
        if sh([os.path.join(ENV, "run_pypeit"), os.path.basename(pfile), "-c"],
              log, cwd=sdir):
            print(f"{tag}: CALIB RUN FAILED, night continues with next setup",
                  flush=True)
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

        # star position -> sky regions, from the LONGEST science exposure.  The
        # finder needs the best signal it can get; file order is meaningless here
        # and on 2023-11-07 it put a 5 s frame first.
        best = max(sci, key=lambda f: exptime.get(f, 0.0))
        out = subprocess.run(
            [os.path.join(ENV, "python"), FINDER,
             os.path.join(RAW, night, best), os.path.join(sdir, "Calibrations")],
            capture_output=True, text=True)
        open(log, "a").write(out.stdout + out.stderr)
        # No source found.  Stop rather than reduce: a degenerate ':0,100:' region
        # leaves PypeIt with no sky pixels, and it will still produce a spec2d --
        # one whose sky model check_skysub.py cannot validate (nan sky-zone
        # residuals, meaningless "sky ok" verdict).  That is a silent bad result,
        # so it needs a human look before anything is spent on it.
        if "NO DETECTION" in out.stdout:
            with open(NEEDS_DECISION, "a") as f:
                f.write(f"{tag}\tno object detected in {best} "
                        f"({exptime.get(best, 0.0):.0f}s)\n")
            print(f"{tag}: *** NO OBJECT DETECTED -- NOT REDUCED, awaiting your "
                  f"confirmation ***", flush=True)
            print(f"{tag}:     checked {best} ({exptime.get(best, 0.0):.0f}s, the "
                  f"longest of {len(sci)} frames)", flush=True)
            print(f"{tag}:     inspect the frame for a source; logged to "
                  f"{os.path.basename(NEEDS_DECISION)}", flush=True)
            continue

        m = re.search(r"user_regions = (\S+)", out.stdout)
        if not m:
            print(f"{tag}: FIND_OBJECT_REGIONS FAILED (no user_regions in its "
                  f"output), setup skipped", flush=True)
            continue
        regions = m.group(1)
        print(f"{tag}: user_regions = {regions}", flush=True)

        txt = open(pfile).read()
        if "[[skysub]]" not in txt:
            txt = txt.replace(
                f"        reid_arxiv = {tmpl}",
                f"        reid_arxiv = {tmpl}\n\n[reduce]\n"
                f"    [[skysub]]\n        user_regions = {regions}", 1)
            open(pfile, "w").write(txt)

        print(f"{tag}: science run ...", flush=True)
        rc = sh([os.path.join(ENV, "run_pypeit"), os.path.basename(pfile)],
                log, cwd=sdir)
        n2d = len(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
        print(f"{tag}: {'DONE' if rc == 0 else f'FAILED rc={rc}'} "
              f"({n2d}/{len(sci)} spec2d)", flush=True)

print("ALL NIGHTS PROCESSED", flush=True)
