#!/usr/bin/env python
"""Drive the RH4 standard-star reduction, with the RL/RH1/RH2 procedure.

RH4 is one night, 2024-12-28, one usable configuration, and **two standards in
it**: feige110 (3 frames, 50-240 s) and feige34 (2 frames, 180 s each).

Per setup:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RH4/<night> -c all   (seed_rh4.py)
  2. ThAr frames retyped arc,tilt and `lamps` pointed at ThArRH4     (seed_rh4.py)
  3. run_pypeit -c   with the template                              (gate_rh4_template.py)
  4. find_object_regions.py PER STAR -> user_regions
  5. run_pypeit                                                     (here)

What differs from run_rh2.py, and why:

**The wavelength settings are a pair, not a template.**  RH4 needs
`reid_arxiv = keck_kcrm_RH4.fits` *and* `lamps = ThArRH4`.  PypeIt's KCWI
catalogue holds four lines above 9550 A, and `full_template` fits to catalogue
wavelengths, so a run that kept the default lamps would fail with the template
in place and look like a bad template.  This driver checks for both and refuses
to run with only one.

**Sky regions are found per STAR, not per setup.**  Both stars are in setup B
because `pypeit_setup` groups on (dispname, decker, binning, cenwave) and the
target is not part of that key.  They were acquired 11 hours apart, so there is
no reason for them to land on the same slices.  The driver measures each, and:

  * if the two regions agree within `--region-tol` percent, it reduces the setup
    in one pass with the shared value;
  * if they do not, it reduces each star in its own pass, with the other star's
    rows commented out, so each gets its own sky definition.

Calibrations are built once and reused across those passes -- PypeIt keys them by
frame, not by which science rows are active.

Retained from run_rh2.py unchanged: the incomplete-config skip (RH4's setup A is
the night's stray DARK plus seven biases, no arc/trace/pixelflat), the
MIN_SCI_EXPTIME floor and longest-exposure frame choice for the finder, the
`NO DETECTION` stop, and the rule that a pre-template `Calibrations/` must be
deleted rather than reused.

Usage:
    python run_rh4.py [--region-tol 3] [--dry-run]
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
RAW = os.path.join(ROOT, "fits", "by_night", "RH4")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
NEEDS_DECISION = os.path.join(BASE, "NEEDS_DECISION.txt")

TEMPLATE = os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH4.fits")
LAMPS = os.path.join(ROOT, "pypeit_test", "ThArRH4")

# Below this an exposure is an acquisition or focus test, not a standard.
MIN_SCI_EXPTIME = 10.0
# RH4's slices span ~985 A, against the 400-900 A pick_seed_slits was given for
# RH1 and RH2.  Span is a property of the grating and has to be restated.
SPAN = ("800", "1200")


def frames(pfile):
    """({frametype: [filename, ...]}, {filename: exptime}, {filename: target}).

    Columns are located by name from the data block's own header line rather
    than by position: a .pypeit written by another PypeIt version can order them
    differently, and reading exptime out of the wrong column silently drops good
    frames.
    """
    cols, types, exp, tgt = None, {}, {}, {}
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
        tgt[fn] = re.sub(r"[^a-z0-9]+", "",
                         (row.get("target", "") or "").strip().lower()) or "std"
    return types, exp, tgt


def sh(cmd, log, cwd=None):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def set_regions(pfile, regions):
    """Write user_regions, replacing any that is already there."""
    txt = open(pfile).read()
    if "[[skysub]]" in txt:
        txt = re.sub(r"user_regions = \S+", f"user_regions = {regions}", txt, count=1)
    else:
        txt = txt.replace(
            f"        lamps = {LAMPS}",
            f"        lamps = {LAMPS}\n\n[reduce]\n"
            f"    [[skysub]]\n        user_regions = {regions}", 1)
    open(pfile, "w").write(txt)


def activate(pfile, keep):
    """Leave only `keep` science rows uncommented; other stars are commented out.

    Calibration rows are never touched, so the same Calibrations/ serves every
    pass.  Idempotent: a row is re-activated if it is in `keep` and commented.
    """
    out = []
    for line in open(pfile).read().split("\n"):
        bare = line[1:].strip() if line.startswith("#") else line
        parts = bare.split("|")
        fn = parts[0].strip()
        is_sci = (fn.endswith(".fits") and len(parts) > 1
                  and parts[1].strip() == "science")
        if not is_sci:
            out.append(line)
        elif fn in keep:
            out.append(bare)
        else:
            out.append("# " + bare if not line.startswith("#") else line)
    open(pfile, "w").write("\n".join(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region-tol", type=float, default=3.0,
                    help="percent of a slice by which two stars' sky regions may "
                         "differ and still share one reduction pass")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(TEMPLATE):
        sys.exit(f"{TEMPLATE} missing -- run build_rh4_template.py")
    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"{LAMPS}_lines.dat missing -- run build_rh4_linelist.py")

    for night in sorted(os.listdir(RAW)):
        ndir = os.path.join(BASE, night)
        if not os.path.isdir(os.path.join(RAW, night)) or not os.path.isdir(ndir):
            continue
        log = os.path.join(ndir, "driver.log")
        rundir = os.path.join(ndir, "pypeit_run")

        for sdir in sorted(glob.glob(os.path.join(rundir, "keck_kcrm_*"))):
            pfile = sorted(glob.glob(os.path.join(sdir, "*.pypeit")))[0]
            tag = f"{night}/{os.path.basename(sdir)}"
            ftypes, exptime, target = frames(pfile)

            sci = ftypes.get("science", [])
            short = [f for f in sci if not (exptime.get(f, 0.0) >= MIN_SCI_EXPTIME)]
            if short:
                print(f"{tag}: ignoring {len(short)} science frame(s) under "
                      f"{MIN_SCI_EXPTIME:.0f}s -- test frames", flush=True)
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
            if TEMPLATE not in txt or LAMPS not in txt:
                print(f"{tag}: .pypeit is not pointed at BOTH keck_kcrm_RH4.fits "
                      f"and ThArRH4 -- run gate_rh4_template.py first, skipped",
                      flush=True)
                continue

            # How many of the 24 slices actually got a wavelength solution.  A
            # setup reduced with slices masked is not the same measurement as one
            # with all of them, and nothing downstream says so out loud.
            wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
            if not wc:
                print(f"{tag}: no WaveCalib -- run gate_rh4_template.py first, "
                      f"skipped", flush=True)
                continue
            chk = subprocess.run([os.path.join(ENV, "python"), SEEDS,
                                  "--span", SPAN[0], SPAN[1], wc[0]],
                                 capture_output=True, text=True)
            open(log, "a").write(chk.stdout + chk.stderr)
            m = re.search(r"(\d+) usable", chk.stdout)
            print(f"{tag}: wavelength solutions on "
                  f"{m.group(1) if m else '?'} slices", flush=True)

            # ---- sky regions, per star ------------------------------------
            stars = sorted({target[f] for f in sci})
            regions, best_of = {}, {}
            for star in stars:
                fs = [f for f in sci if target[f] == star]
                best = max(fs, key=lambda f: exptime.get(f, 0.0))
                best_of[star] = best
                out = subprocess.run(
                    [os.path.join(ENV, "python"), FINDER,
                     os.path.join(RAW, night, best),
                     os.path.join(sdir, "Calibrations")],
                    capture_output=True, text=True)
                open(log, "a").write(out.stdout + out.stderr)
                if "NO DETECTION" in out.stdout:
                    with open(NEEDS_DECISION, "a") as f:
                        f.write(f"{tag}\t{star}: no object detected in {best} "
                                f"({exptime.get(best, 0.0):.0f}s)\n")
                    print(f"{tag}: *** NO OBJECT DETECTED for {star} -- NOT "
                          f"REDUCED ***", flush=True)
                    continue
                mm = re.search(r"user_regions = (\S+)", out.stdout)
                if not mm:
                    print(f"{tag}: FIND_OBJECT_REGIONS FAILED for {star}",
                          flush=True)
                    continue
                regions[star] = mm.group(1)
                print(f"{tag}: {star} ({exptime.get(best,0):.0f}s) "
                      f"user_regions = {regions[star]}", flush=True)

            if not regions:
                print(f"{tag}: no usable sky regions, setup skipped", flush=True)
                continue

            def bounds(r):
                a, b = r.split(",")
                return float(a.lstrip(":") or 0), float(b.rstrip(":") or 100)

            vals = [bounds(r) for r in regions.values()]
            spread = max(max(abs(v[0] - w[0]), abs(v[1] - w[1]))
                         for v in vals for w in vals)
            shared = spread <= args.region_tol
            print(f"{tag}: sky regions across {len(regions)} star(s) differ by "
                  f"{spread:.1f}% -> {'ONE pass' if shared else 'ONE PASS PER STAR'}",
                  flush=True)

            if args.dry_run:
                continue

            passes = ([(None, list(regions.values())[0], sci)] if shared else
                      [(s, regions[s], [f for f in sci if target[f] == s])
                       for s in sorted(regions)])

            for star, reg, keep in passes:
                label = star or "all stars"
                activate(pfile, set(keep))
                set_regions(pfile, reg)
                print(f"{tag}: science run [{label}], {len(keep)} frames, "
                      f"user_regions = {reg} ...", flush=True)
                rc = sh([os.path.join(ENV, "run_pypeit"),
                         os.path.basename(pfile)], log, cwd=sdir)
                n2d = len(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
                print(f"{tag}: [{label}] {'DONE' if rc == 0 else f'FAILED rc={rc}'} "
                      f"({n2d} spec2d so far)", flush=True)

            # Leave every science row active, so the file describes the night
            # rather than the last pass over it.
            activate(pfile, set(sci))

    print("ALL NIGHTS PROCESSED", flush=True)


if __name__ == "__main__":
    main()
