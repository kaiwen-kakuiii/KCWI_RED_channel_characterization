#!/usr/bin/env python
"""Drive the RH3 standard-star reduction, with the RL/RH1/RH2/RH4 procedure.

RH3 is four nights, one usable configuration each (setup B; setup A is always
that night's stray DARK plus biases), Medium slicer and 2x2 binning throughout:

    2023-11-12 B   cenwave 8600   feige110   2 frames  (180, 360 s)
    2024-10-01 B   cenwave 8400   g191b2b    4 frames  (1, 100, 100, 100 s)
    2024-10-29 B   cenwave 8600   feige110   2 frames  (60, 300 s)
    2024-11-06 B   cenwave 8600   feige110   2 frames  (300, 300 s)

Per setup:

  1. pypeit_setup -s keck_kcrm -r fits/by_night/RH3/<night> -c all
  2. holy-grail seeds, default FeAr lamps                    (seed_rh3.py)
  3. build one template per CENWAVE                          (build_rh3_template.py)
  4. run_pypeit -c against it                                (gate_rh3_template.py)
  5. find_object_regions.py -> user_regions
  6. run_pypeit                                              (here)

What differs from run_rh4.py, and why:

**The template and the ThAr line list are one unit, as on RH4.**  RH3's arc is
not FeAr: `seed_rh3.py` retypes the ThAr frames as the arc and comments the FeAr
frames out, because FeAr's 23 in-band catalogue lines are too few for holy-grail
to identify a solution from scratch.  With a ThAr arc and a ThAr-derived
template, the catalogue must be ThAr as well -- measured on 2023-11-12, the
default FeI/ArI/ArII list fits 9-12 lines per slit against ThArRH3's 58-78,
because only the argon lines of a ThAr spectrum appear in it.  That still
reported 24/24 at rms 0.029 px, which is the "few points, small residual" trap,
not a good solution.  So this driver **refuses to start unless BOTH the cenwave-
matched template and `lamps = ThArRH3` are set**, the same refusal run_rh4.py
makes for the same reason.

**The template is chosen per setup, by cenwave.**  Three nights sit at 8600 and
one at 8400, and a template 200 A off is what broke the shipped
`keck_kcrm_RH3.fits` in the first place.  Handing a setup the wrong one would
reproduce that failure exactly, so the cenwave in the .pypeit is read and matched
rather than assumed from the night.

**One star per configuration**, unlike RH4, so the per-star reduction passes
never fire.  The machinery is kept because a driver that silently coadds two
different stars is the failure RH4 documents, and nothing guarantees the next
RH3 night will be single-star.

Retained unchanged from run_rh4.py: the incomplete-config skip, the
MIN_SCI_EXPTIME floor and longest-exposure frame choice for the finder (which
drops 2024-10-01's 1 s acquisition frame), the `NO DETECTION` stop, and the rule
that a pre-template `Calibrations/` must be deleted rather than reused.

Usage:
    python run_rh3.py [--region-tol 3] [--dry-run]
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
RAW = os.path.join(ROOT, "fits", "by_night", "RH3")
FINDER = os.path.join(ROOT, "pypeit_test", "find_object_regions.py")
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
NEEDS_DECISION = os.path.join(BASE, "NEEDS_DECISION.txt")

TEMPLATE = {8600: os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH3_8600.fits"),
            8400: os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH3_8400.fits")}
LAMPS = os.path.join(ROOT, "pypeit_test", "ThArRH3")

# Below this an exposure is an acquisition or focus test, not a standard.
# 2024-10-01's 1 s g191b2b frame is the one this drops.
MIN_SCI_EXPTIME = 10.0

# RH3 slices span ~793 A, inside pick_seed_slits' 400-900 A default, so no
# --span override is needed.  Stated because RH4 lost ten minutes to a window
# set from a different grating.
SPAN = None


def cenwave_of(pfile):
    m = re.search(r"cenwave:\s*([0-9.]+)", open(pfile).read())
    return None if not m else int(round(float(m.group(1)) / 50.0) * 50)


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

    The whole header is rebuilt rather than spliced into, because splicing broke:
    the previous version inserted the `[reduce]` block immediately after the
    `reid_arxiv` line, which was the last line of the wavelengths block until
    `lamps` was added below it.  After that, inserting there cut the block in two
    and left `lamps` sitting under `[[skysub]]`, so every science run died with

        ValueError: ['lamps'] not recognized key(s) for SkySubPar.

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
    lamps = value("lamps")
    if reid is None:
        raise RuntimeError(f"{pfile}: no reid_arxiv -- run gate_rh3_template.py")
    if lamps is None:
        raise RuntimeError(f"{pfile}: no lamps -- the arc is ThAr, so the "
                           f"catalogue must be ThArRH3; run gate_rh3_template.py")

    head = ("# Auto-generated PypeIt input file using PypeIt version: 2.0.1\n"
            "\n# User-defined execution parameters\n"
            "[rdx]\n    spectrograph = keck_kcrm\n\n"
            "[calibrations]\n    [[wavelengths]]\n"
            f"        method = {method}\n"
            f"        reid_arxiv = {reid}\n"
            f"        lamps = {lamps}\n\n"
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
    ap.add_argument("--only", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # Both halves, checked before anything runs.  A setup pointed at the template
    # but left on the default lamps fails in a way that looks like a bad
    # template -- it reports 24/24 at a low RMS through ~11 fitted lines.
    missing = [f for f in (LAMPS + "_lines.dat",) if not os.path.exists(f)]
    missing += [t for t in TEMPLATE.values() if not os.path.exists(t)]
    if len(missing) == len(TEMPLATE) + 1:
        sys.exit("neither template nor line list exists -- run "
                 "build_rh3_template.py and build_rh4_linelist.py --label RH3")
    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"{LAMPS}_lines.dat missing -- the arc is ThAr, so the "
                 f"catalogue must be ThAr")

    for night in sorted(os.listdir(RAW)):
        if args.only and args.only != night:
            continue
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

            cw = cenwave_of(pfile)
            if cw not in TEMPLATE:
                print(f"{tag}: cenwave {cw} has no RH3 template, skipped", flush=True)
                continue
            txt = open(pfile).read()
            if TEMPLATE[cw] not in txt:
                print(f"{tag}: .pypeit is not pointed at "
                      f"{os.path.basename(TEMPLATE[cw])} -- run "
                      f"gate_rh3_template.py --night {night} --cenwave {cw} "
                      f"first, skipped", flush=True)
                continue
            if LAMPS not in txt:
                print(f"{tag}: .pypeit is not pointed at ThArRH3.  The arc frames "
                      f"are ThAr, so PypeIt's default FeI/ArI/ArII catalogue "
                      f"would fit ~11 lines per slit instead of ~65.  Run "
                      f"gate_rh3_template.py first.  Skipped.", flush=True)
                continue

            # How many of the 24 slices actually got a wavelength solution.  A
            # setup reduced with slices masked is not the same measurement as one
            # with all of them, and nothing downstream says so out loud.
            wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
            if not wc:
                print(f"{tag}: no WaveCalib -- run gate_rh3_template.py first, "
                      f"skipped", flush=True)
                continue
            cmd = [os.path.join(ENV, "python"), SEEDS]
            if SPAN:
                cmd += ["--span", SPAN[0], SPAN[1]]
            chk = subprocess.run(cmd + [wc[0]], capture_output=True, text=True)
            open(log, "a").write(chk.stdout + chk.stderr)
            m = re.search(r"(\d+) usable", chk.stdout)
            print(f"{tag}: cenwave {cw}, wavelength solutions on "
                  f"{m.group(1) if m else '?'} slices", flush=True)

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

    print("ALL NIGHTS PROCESSED", flush=True)


if __name__ == "__main__":
    main()
