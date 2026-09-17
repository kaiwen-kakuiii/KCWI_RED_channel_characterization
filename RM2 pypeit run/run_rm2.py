#!/usr/bin/env python
"""Drive the RM2 standard-star reduction, with the RL/RH1/RH2/RH3/RH4/RM1 procedure.

RM2 is nine nights carrying nine science configurations at five central
wavelengths, in two slicers and TWO STANDARD STARS:

    2024-05-09 A   cenwave 7750   Small  1x1   3 frames  (120 s each)      feige34
    2024-06-10 A   cenwave 8850   Small  1x1   6 frames  (10 s each)       feige110
    2024-12-24 A   cenwave 8850   Small  1x1   4 frames  (10, 20, 20, 20)  feige 34
    2024-01-04 B   cenwave 8900   Medium 2x2   3 frames  (120 s each)      feige34
    2024-04-30 B   cenwave 8900   Medium 2x2   3 frames  (150 s each)      feige34
    2024-06-11 B   cenwave 8950   Medium 2x2   4 frames  (25, 40, 40, 40)  feige34
    2025-01-01 B   cenwave 8950   Medium 2x2   4 frames  (66, 68, 35, 35)  feige34
    2025-01-02 B   cenwave 8950   Medium 2x2   3 frames  (35 s each)       feige34
    2023-09-23 B   cenwave 9850   Medium 2x2   6 frames  (3, 30, 60x4)     feige110

Per setup, unchanged from run_rm2.py:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RM2/<night> -c all   (setup_rm2.sh)
  2. gate the registered template                      (gate_rm2_template.py)
     -- building one first where the shipped template fails
     (build_rm2_template.py)
  3. find_object_regions.py -> user_regions
  4. run_pypeit                                        (here)

## What differs from run_rm2.py, and why

**One setup per night, and it is not always B.**  Every night here carries
exactly one science configuration; the other setup, where there is one, is the
night's stray FPCam/dark/bias group at cenwave 0 or 24009 with zero science
frames (RM1 section 1, RH3's "setup A is always that night's stray DARK plus
biases").  On 2024-05-09, 2024-06-10 and 2024-12-24 the SCIENCE setup is
lettered A because those nights have no stray group at all.  Filter on the
science-frame count, never on the letter -- the letter is not stable.

**Three nights at one central wavelength.**  8950 has 2024-06-11, 2025-01-01 and
2025-01-02; 8900 and 8850 have two nights each.  That is the strongest
reproducibility test in this project so far, and what settles polynomial order
(order_study_rm2.py) instead of inheriting RM1's 5.

**Two stars, and at 8850 they differ BETWEEN the two nights** -- feige110 on
2024-06-10, feige 34 (with a space, as the header writes it) on 2024-12-24.  No
single setup mixes stars, so the per-star reduction passes never fire, but the
8850 pair is a cross-STAR comparison and not only a cross-night one.

**A slice spans ~1994 A, not RM1's ~1453.**  `pick_seed_slits` is called with
--span 1700 2300; on its 400-900 default, or on RM1's 1200-1700, every RM2 slice
is rejected as nonsense and a perfect setup reports 0 usable.  Same trap, third
grating, different numbers.

**The SHIFTED tolerance is 150/180 A, measured not inherited.**  On 2023-09-23 B
-- gated 24/24, dispersions agreeing to 0.4% -- the 24 blue ends span 206.7 A
with a worst deviation of 128.4 A from the median, in a clean odd/even comb.
RM1's Large 100 A would flag that correct solution.  See gate_rm2_template.py.

**Line counts are thin and must be watched.**  RM1's FeAr arc fitted 70-80 lines
per slit; RM2's 9850 setup fits 10-14.  The catalogue is the same PypeIt default
FeI/ArI/ArII list, and RM2 sits 1000-3000 A redder, where it runs out.  A fit
through 12 lines returns a small rms (0.026 px here) whether or not it is right
-- the "few points, small residual" trap of RH1_PROCEDURE.md -- so the check
that matters is dispersion agreement across the 24 slices and, at 8850/8900/8950,
agreement between independent nights.

Retained unchanged from run_rm2.py: the 10 s exposure floor that yields at two
frames, the incomplete-config skip, the longest-exposure frame choice for the
finder, the `NO DETECTION` stop, the per-star reduction passes, and the rule that
a pre-gate `Calibrations/` must be deleted rather than reused.

Usage:
    python run_rm2.py [--region-tol 3] [--only 2024-04-30] [--dry-run]
"""

import argparse
import glob
import os
import re
import subprocess
import sys

from astropy.io import fits

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
RAW = os.path.join(ROOT, "fits", "by_night", "RM2")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
NEEDS_DECISION = os.path.join(BASE, "NEEDS_DECISION.txt")

# One source of truth for which template serves which central wavelength: the
# gate decides it, the driver only checks the .pypeit actually says so.  Two
# copies of this mapping would eventually disagree, and the failure mode --
# reducing against a template 200 A off -- is silent.
sys.path.insert(0, os.path.join(ROOT, "pypeit_test"))
from gate_rm2_template import (REGISTRY, SHIPPED, reid_value,   # noqa: E402
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

# Physically possible span of ONE RM2 slice, in A: 2064 binned px x 0.9659,
# measured 1958-1973 on 2023-09-23 B.
SPAN = ("1700", "2300")

# Both windows must be passed to pick_seed_slits, not just the span.  Left on its
# 50 A default, `shift_tol` flags RM1's genuine 84 A slice-to-slice spread as
# SHIFTED and this driver reported "wavelength solutions on 20 slices" for a
# setup the gate had passed at 24/24.  Nothing downstream depended on that
# number, which is exactly why it would have been believed.
# Picked per SLICER, and RM2's slicers are Medium and Small -- there is no Large
# config here, so RM1's values cannot be carried over by decker name.  Measured
# 128.4 A from the median on 2023-09-23 B (Medium), a setup gated at 24/24.
def shift_tol_of(pfile):
    m = re.search(r"decker:\s*(\S+)", open(pfile).read())
    return str(shift_tol_for(m.group(1) if m else "Medium"))


def cenwave_of(pfile):
    """Rounded to 10 A: RM2's 8850/8900/8950 settings sit 50 A apart and a
    coarser round would merge three different grating angles into one."""
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
        raise RuntimeError(f"{pfile}: no reid_arxiv -- run gate_rm2_template.py")

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
                      f"gate_rm2_template.py --night {night} --setup "
                      f"{os.path.basename(sdir)} first, skipped", flush=True)
                continue

            # How many of the 24 slices actually got a wavelength solution.  A
            # setup reduced with slices masked is not the same measurement as one
            # with all of them, and nothing downstream says so out loud.
            wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
            if not wc:
                print(f"{tag}: no WaveCalib -- run gate_rm2_template.py first, "
                      f"skipped", flush=True)
                continue
            chk = subprocess.run([os.path.join(ENV, "python"), SEEDS,
                                  "--span", SPAN[0], SPAN[1],
                                  "--shift-tol", shift_tol_of(pfile),
                                  "--nline-min", "8", wc[0]],
                                 capture_output=True, text=True)
            open(log, "a").write(chk.stdout + chk.stderr)
            m = re.search(r"(\d+) usable", chk.stdout)
            print(f"{tag}: cenwave {cw}, template "
                  f"{os.path.basename(REGISTRY[cw])}, wavelength solutions on "
                  f"{m.group(1) if m else '?'} of 24 slices usable", flush=True)

            # ---- sky regions, per star AND PER POINTING -------------------
            # RM1's driver derived ONE region per star per configuration, from
            # that star's longest exposure, and reduced every frame against it.
            # RM2 2023-09-23 B is where that breaks: its five frames sit at
            # (0.0, 0.0) x3 and (-0.7, +1.0) x2 -- a 1.2 arcsec nod -- and
            # find_object_regions, run on each frame separately, does not agree
            # across it:
            #
            #     frames 1-3   brightest slice 12   object 44-57%   -> :38,65:
            #     frames 4-5   brightest slice 11   object 38-51%   -> :31,59:
            #
            # Reduced against the pointing-A window, frames 4-5 have the star's
            # blue wing sitting inside the region declared SKY, so the sky
            # b-spline fits part of the star and subtracts it.  Measured on those
            # frames: star flux DOWN 17% and sky flux UP 12% across the nod,
            # anticorrelated, which is the signature.  It cost 9850 most of its
            # throughput deficit and looked exactly like cloud.
            #
            # `user_regions` is one parameter for a whole run_pypeit call, so the
            # only way to honour two windows is two passes.  The machinery for
            # that already existed for two STARS; the key is now (star, pointing)
            # and everything downstream is unchanged.
            def pointing_of(fn):
                h = fits.getheader(os.path.join(RAW, night, fn))
                return (round(float(h.get("RAOFF", 0.0) or 0.0), 1),
                        round(float(h.get("DECOFF", 0.0) or 0.0), 1))

            pt = {f: pointing_of(f) for f in sci}
            keys = sorted({(target[f], pt[f]) for f in sci})
            if len({k[1] for k in keys}) > 1:
                print(f"{tag}: {len({k[1] for k in keys})} pointings "
                      f"{sorted({k[1] for k in keys})} -- sky regions measured "
                      f"per pointing", flush=True)

            regions = {}
            for key in keys:
                star, point = key
                fs = [f for f in sci if (target[f], pt[f]) == key]
                best = max(fs, key=lambda f: exptime.get(f, 0.0))
                out = subprocess.run(
                    [os.path.join(ENV, "python"), FINDER,
                     os.path.join(RAW, night, best),
                     os.path.join(sdir, "Calibrations")],
                    capture_output=True, text=True)
                open(log, "a").write(out.stdout + out.stderr)
                if "NO DETECTION" in out.stdout:
                    with open(NEEDS_DECISION, "a") as f:
                        f.write(f"{tag}\t{star} at {point}: no object detected "
                                f"in {best} ({exptime.get(best, 0.0):.0f}s)\n")
                    print(f"{tag}: *** NO OBJECT DETECTED for {star} at {point} "
                          f"-- NOT REDUCED ***", flush=True)
                    continue
                mm = re.search(r"user_regions = (\S+)", out.stdout)
                if not mm:
                    print(f"{tag}: FIND_OBJECT_REGIONS FAILED for {star} at "
                          f"{point}", flush=True)
                    continue
                regions[key] = mm.group(1)
                print(f"{tag}: {star} at {point} ({exptime.get(best,0):.0f}s) "
                      f"user_regions = {regions[key]}", flush=True)

            if not regions:
                print(f"{tag}: no usable sky regions, setup skipped", flush=True)
                continue

            def bounds(r):
                a, b = r.split(",")
                return float(a.lstrip(":") or 0), float(b.rstrip(":") or 100)

            # Cluster the (star, pointing) groups whose windows agree to within
            # --region-tol, so a nod too small to move the star does not buy an
            # extra reduction pass and half the frames per cube.
            clusters = []          # [[key, ...], ...], first key carries the window
            for key in sorted(regions):
                v = bounds(regions[key])
                for cl in clusters:
                    w = bounds(regions[cl[0]])
                    if max(abs(v[0] - w[0]), abs(v[1] - w[1])) <= args.region_tol:
                        cl.append(key)
                        break
                else:
                    clusters.append([key])

            spread = max((max(abs(bounds(regions[a])[0] - bounds(regions[b])[0]),
                              abs(bounds(regions[a])[1] - bounds(regions[b])[1]))
                          for a in regions for b in regions), default=0.0)
            print(f"{tag}: sky regions across {len(regions)} (star, pointing) "
                  f"group(s) differ by {spread:.1f}% -> "
                  f"{len(clusters)} reduction pass(es)", flush=True)

            if args.dry_run:
                continue

            passes = []
            for cl in clusters:
                keep = [f for f in sci if (target[f], pt[f]) in set(cl)]
                label = ", ".join(f"{k[0]}@{k[1]}" for k in cl)
                passes.append((label, regions[cl[0]], keep))

            for label, reg, keep in passes:
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
