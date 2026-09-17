#!/usr/bin/env python
"""Drive the RM1 standard-star reduction, with the RL/RH1/RH2/RH3/RH4 procedure.

RM1 is six nights carrying **nine** science configurations, every one of them
feige34:

    2024-04-01 A   cenwave 6130   Small 1x1   2 frames  ( 25, 100 s)
    2024-05-02 B   cenwave 6200   Large 2x2   6 frames  (5, 15, 45, 60, 120, 120 s)
    2024-03-15 B   cenwave 6300   Large 2x2   7 frames  (6, 15, 45 x5)
    2023-12-09 B   cenwave 6480   Large 2x2   2 frames  (  7,  14 s)
    2023-12-10 C   cenwave 6630   Large 2x2   4 frames  (15 s each)
    2023-12-11 C   cenwave 6630   Large 2x2   4 frames  (15 s each)
    2024-04-01 B   cenwave 7390   Small 1x1   2 frames  (200 s each)
    2023-12-10 B   cenwave 7510   Large 2x2   3 frames  (15 s each)
    2023-12-11 B   cenwave 7510   Large 2x2   3 frames  (15 s each)

Per setup:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RM1/<night> -c all
  2. gate the registered template                      (gate_rm1_template.py)
     -- building one first where the shipped template fails
     (build_rm1_template.py)
  3. find_object_regions.py -> user_regions
  4. run_pypeit                                        (here)

Holy-grail never enters, unlike RH1/RH2/RH3: PypeIt ships an RM1 template that
genuinely works at three of the seven central wavelengths, and the four it does
not reach are reached by bootstrapping from its own good slices.

What differs from run_rh3.py, and why:

**A night is not a configuration here.**  2023-12-10 and 2023-12-11 each carry
TWO science setups -- B at cenwave 7510 and C at 6630 -- so everything keys on
the setup directory and nothing keys on the date.  run_rh3.py already looped over
setups; RM1 is the first run where that loop does more than one useful pass.

**No custom line list.**  RH3 and RH4 had to retype the ThAr frames and build
their own catalogues because the default FeI/ArI/ArII list was too thin in band.
RM1's FeAr arc is well served: measured on 2023-12-10 C, the default catalogue
fitted 70-80 lines per slit at rms 0.200 px.  So this driver checks only that
the .pypeit points at the right TEMPLATE, and deliberately does not demand a
`lamps` line the way run_rh3.py must.

**The exposure-time floor is 3 s, not RH3's 10 s.**  That rule exists to drop
acquisition and focus frames -- RH3's 1 s g191b2b frame.  RM1's shortest science
frames are 5, 6 and 7 s, and they are real exposures of a bright standard
(feige34, V = 11.2) at medium dispersion.  Applying RH3's 10 s floor here would
drop the 7 s frame of 2023-12-09 B, leaving that configuration with ONE frame --
and a single-frame cube is the unsolved failure THROUGHPUT_PROCEDURE.md
documents, where `combine = False` writes a broken WCS and inflates the counts by
~1e21.  A floor imported from another grating would have destroyed a whole
central wavelength.

**pick_seed_slits needs --span.**  An RM1 slice spans ~1453 A, far outside the
400-900 A default set from RH1 and RH2, and with the default every slice is
rejected as nonsense.  RH3's notes warn about exactly this after RH4 lost ten
minutes to a window set from a different grating; RM1 is where it would have
been silent rather than obvious, because the reduction still runs.

Retained unchanged from run_rh3.py: the incomplete-config skip, the longest-
exposure frame choice for the finder, the `NO DETECTION` stop, the per-star
reduction passes (RM1 is single-star, so they never fire, but a driver that
silently coadds two stars is the failure RH4 documents), and the rule that a
pre-gate `Calibrations/` must be deleted rather than reused.

Usage:
    python run_rm1.py [--region-tol 3] [--only 2024-03-15] [--dry-run]
"""
import argparse
import glob
import os
import re
import subprocess
import sys

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
RAW = os.path.join(ROOT, "fits", "by_night", "RM1")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
NEEDS_DECISION = os.path.join(BASE, "NEEDS_DECISION.txt")

# One source of truth for which template serves which central wavelength: the
# gate decides it, the driver only checks the .pypeit actually says so.  Two
# copies of this mapping would eventually disagree, and the failure mode --
# reducing against a template 200 A off -- is silent.
sys.path.insert(0, os.path.join(ROOT, "pypeit_test"))
from gate_rm1_template import (REGISTRY, SHIPPED, reid_value,   # noqa: E402
                               shift_tol_for)

# Below this an exposure is an acquisition or focus test, not a standard.  The
# project rule is 10 s: short frames sit at the START of a sequence, before
# guiding and focus settle, so exposure time and frame quality are correlated --
# measured on 2024-03-15 B, whose 6 s and 15 s frames were taken at GUIDFWHM
# 2.21" and 1.77" against 1.18-1.34" for the five 45 s frames that followed, and
# came out ~12% low.
#
# But the floor YIELDS when honouring it would leave a configuration with fewer
# than MIN_FRAMES_AFTER_CUT frames.  A single-frame cube writes a broken WCS and
# inflates counts by ~1e21 (THROUGHPUT_PROCEDURE.md), so enforcing the rule there
# destroys the measurement it was meant to protect.  RM1 2023-12-09 B is exactly
# that case: 7 s + 14 s, and it is the only coverage at cenwave 6480.
MIN_SCI_EXPTIME = 10.0
MIN_FRAMES_AFTER_CUT = 2


def apply_exptime_floor(sci, exptime, tag):
    """Drop sub-floor frames, unless dropping them would starve the config."""
    short = [f for f in sci if exptime.get(f, 0.0) < MIN_SCI_EXPTIME]
    keep = [f for f in sci if f not in short]
    if not short:
        return sci
    if len(keep) >= MIN_FRAMES_AFTER_CUT:
        print(f"{tag}: dropping {len(short)} frame(s) under "
              f"{MIN_SCI_EXPTIME:.0f}s, {len(keep)} remain", flush=True)
        return keep
    print(f"{tag}: KEEPING {len(short)} frame(s) under {MIN_SCI_EXPTIME:.0f}s -- "
          f"dropping them would leave {len(keep)}, under the {MIN_FRAMES_AFTER_CUT} "
          f"needed for a valid cube", flush=True)
    return sci

# Physically possible span of ONE RM1 slice, in A.  Measured ~1453.
SPAN = ("1200", "1700")

# Both windows must be passed to pick_seed_slits, not just the span.  Left on its
# 50 A default, `shift_tol` flags RM1's genuine 84 A slice-to-slice spread as
# SHIFTED and this driver reported "wavelength solutions on 20 slices" for a
# setup the gate had passed at 24/24.  Nothing downstream depended on that
# number, which is exactly why it would have been believed.
# Picked per SLICER: RM1's genuine spread is 79-84 A on Large/2x2 but 106 A on
# Small/1x1, so one number for the grating flags a correct Small solution.
def shift_tol_of(pfile):
    m = re.search(r"decker:\s*(\S+)", open(pfile).read())
    return str(shift_tol_for(m.group(1) if m else "Large"))


def cenwave_of(pfile):
    """Rounded to 10 A, not RH3's 50: RM1's settings sit 100-130 A apart."""
    m = re.search(r"cenwave:\s*([0-9.]+)", open(pfile).read())
    return None if not m else int(round(float(m.group(1)) / 10.0) * 10)


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
    """Rewrite the parameter header canonically, with user_regions set.

    The whole header is rebuilt rather than spliced into, because splicing broke
    on RH3: inserting `[reduce]` immediately after the `reid_arxiv` line cut the
    wavelengths block in two once a `lamps` line sat below it, and every science
    run died with `ValueError: ['lamps'] not recognized key(s) for SkySubPar`.
    Rebuilding is also self-repairing: a file already damaged that way is read
    for its VALUES, wherever they sit, and written back out in the right shape.
    """
    txt = open(pfile).read()

    # Everything from the setup block onward is data and is never touched.
    m = re.search(r"\n(# Setup\n)?setup read\n", txt)
    if not m:
        raise RuntimeError(f"{pfile}: no 'setup read' block found")
    body = txt[m.start():]

    def value(key):
        mm = re.search(rf"^\s*{key}\s*=\s*(\S+)\s*$", txt[:m.start()], re.M)
        return None if mm is None else mm.group(1)

    method = value("method") or "full_template"
    reid = value("reid_arxiv")
    exclude = value("exclude_regions")
    if reid is None:
        raise RuntimeError(f"{pfile}: no reid_arxiv -- run gate_rm1_template.py")

    # `exclude_regions` is carried through for the same reason the gate carries
    # it: it is a hand-made slit-tracing decision (2024-04-01 B), and dropping it
    # here would reduce the science against a different slit set from the one
    # that was gated.
    head = ("# Auto-generated PypeIt input file using PypeIt version: 2.0.1\n"
            "\n# User-defined execution parameters\n"
            "[rdx]\n    spectrograph = keck_kcrm\n\n"
            "[calibrations]\n")
    if exclude:
        head += f"    [[slitedges]]\n        exclude_regions = {exclude}\n"
    head += ("    [[wavelengths]]\n"
             f"        method = {method}\n"
             f"        reid_arxiv = {reid}\n\n"
             "[reduce]\n    [[skysub]]\n"
             f"        user_regions = {regions}\n")
    open(pfile, "w").write(head + body)


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
    ap.add_argument("--only", default=None,
                    help="substring a setup path must contain, e.g. a night or "
                         "'2023-12-10/pypeit_run/keck_kcrm_C' for one config")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for night in sorted(os.listdir(RAW)):
        ndir = os.path.join(BASE, night)
        if not os.path.isdir(os.path.join(RAW, night)) or not os.path.isdir(ndir):
            continue
        log = os.path.join(ndir, "driver.log")
        rundir = os.path.join(ndir, "pypeit_run")

        for sdir in sorted(glob.glob(os.path.join(rundir, "keck_kcrm_*"))):
            if args.only and args.only not in sdir:
                continue
            pfile = sorted(glob.glob(os.path.join(sdir, "*.pypeit")))[0]
            tag = f"{night}/{os.path.basename(sdir)}"
            ftypes, exptime, target = frames(pfile)

            sci = apply_exptime_floor(ftypes.get("science", []), exptime, tag)
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

            cw = cenwave_of(pfile)
            if cw not in REGISTRY:
                print(f"{tag}: cenwave {cw} is not in the template registry, "
                      f"skipped", flush=True)
                continue
            want = reid_value(REGISTRY[cw])
            txt = open(pfile).read()
            if want not in txt:
                print(f"{tag}: .pypeit is not pointed at "
                      f"{os.path.basename(REGISTRY[cw])} -- run "
                      f"gate_rm1_template.py --night {night} --setup "
                      f"{os.path.basename(sdir)} first, skipped", flush=True)
                continue

            # How many of the 24 slices actually got a wavelength solution.  A
            # setup reduced with slices masked is not the same measurement as one
            # with all of them, and nothing downstream says so out loud.
            wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
            if not wc:
                print(f"{tag}: no WaveCalib -- run gate_rm1_template.py first, "
                      f"skipped", flush=True)
                continue
            chk = subprocess.run([os.path.join(ENV, "python"), SEEDS,
                                  "--span", SPAN[0], SPAN[1],
                                  "--shift-tol", shift_tol_of(pfile), wc[0]],
                                 capture_output=True, text=True)
            open(log, "a").write(chk.stdout + chk.stderr)
            m = re.search(r"(\d+) usable", chk.stdout)
            print(f"{tag}: cenwave {cw}, template "
                  f"{os.path.basename(REGISTRY[cw])}, wavelength solutions on "
                  f"{m.group(1) if m else '?'} of 24 slices usable", flush=True)

            # ---- sky regions, per star ------------------------------------
            stars = sorted({target[f] for f in sci})
            regions = {}
            for star in stars:
                fs = [f for f in sci if target[f] == star]
                best = max(fs, key=lambda f: exptime.get(f, 0.0))
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

    print("ALL SETUPS PROCESSED", flush=True)


if __name__ == "__main__":
    main()
