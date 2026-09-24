#!/usr/bin/env python
"""Bootstrap wavelength seeds for RH3, which needs BOTH a template and a lamp swap.

RH3 is the first grating in this project to fail for two reasons at once, and
untangling them took a wrong turn worth recording.

    RH1, RH2   no template shipped at all; FeAr catalogue fine
    RH4        no LINE CATALOGUE in band (6 lines) -- no template could help
    RH3        a template IS shipped but is at the wrong central wavelength,
               AND the FeAr catalogue, while adequate to FIT, is too thin for
               holy-grail to IDENTIFY from scratch

## 1. The shipped template is at the wrong cenwave

`keck_kcwi.py::config_specific_par` maps dispname RH3 -> `keck_kcrm_RH3.fits`, and
PypeIt lists RH3 as a supported setup.  That template is real; it just does not
cover this data:

    shipped keck_kcrm_RH3.fits    8599.4 - 9389.0 A
    cenwave 8600 nights             ~8204 - 8996 A     ~50% overlap
    cenwave 8400 night              ~8004 - 8796 A     ~25% overlap

`full_template` reports the consequence once you read past the failure line:

    Shift = 271.45; cc = 0.1466      slit 1
    Shift = 479.57; cc = 0.2003      slit 2
    Shift = 898.66; cc = 0.1515      slit 7      <- wanders 600 px between slits
    ... Not enough useful IDs                       on all 24

Adjacent KCRM slices differ by a few pixels, not hundreds.  Gate result: 1 usable
slice of 24 at cenwave 8600 and 10 at 8400, every one of them at rms 0.000
through too few points -- 0 usable seeds either way.

## 2. FeAr is adequate to fit but not to identify -- and that distinction cost a run

RH4's rule is to check catalogue coverage before blaming the data.  Done here, it
looks reassuring, and reading it as sufficient was the wrong turn:

    band 8204-8996 A     FeI 6   ArI 14   ArII 3   = 23 lines, median spacing 25 A
    detected, middle slit, on PypeIt's own processed frames:
        FeAr   32 lines >5 sigma, 14 >20 sigma, 7 >100 sigma
        ThAr   78 lines >5 sigma, 35 >20 sigma, 14 >100 sigma

23 catalogue against 32 detected is not the RH4 pathology (6 against 135), and 23
lines *are* enough to constrain the order-4 polynomial `full_template` fits.  So
holy-grail was run on FeAr with the default lamps.  It solved **0 of 24**:

    slit  spat     rms  nlines     wmin      wmax     span     disp
       1   147   0.000       6  -85773.4   58036.1 143809.5  46.2215
       9   747   0.000       5  -25161.5   24713.9  49875.5 -17.4558
      20  1765   0.000       6  -18252.2   16750.5  35002.7 -11.7349

Five to nine lines identified per slit, negative wavelengths, 46 A/px.  Worse than
RH1 (1/24), RH2 (7/24) and RH4 (2/24).

**Fitting and identifying are different requirements, and a catalogue can be
sufficient for one and not the other.**  `full_template` is handed the answer --
the template says which catalogue line each feature is -- and then fits, so 23
lines suffice.  Holy-grail has no template and must recognise the *pattern* of
line spacings by brute force; 23 lines over 793 A, median spacing 25 A and a
112 A gap, do not form patterns distinctive enough to lock onto.  A line count
that clears RH4's bar therefore says nothing about whether holy-grail can work.

## Fix: seed off ThAr, exactly as RH4 does

KCWI takes both lamps and types FeAr `arc`, ThAr `tilt` (`check_frame_type`:
"PypeIt is only setup to wavelength calibrate using the FeAr lamp").  This script
retypes the three ThAr frames `arc,tilt`, comments out the three FeAr frames, and
points `lamps` at a band-limited ThAr list.

`build_rh4_linelist.py --label RH3 --wmin 7900 --wmax 9100 --amp-min 300` writes
`ThArRH3_lines.dat`:

    184 lines, 7901.8-9097.3 A, median spacing 4.27 A, max gap 37 A
    Th I 158,  Th II 16,  Ar I 10

Chosen against what is actually on the frame: 78 detected lines over 793 A is one
per 10.2 A, so a 4.27 A catalogue offers about 2.4 candidates per detected line --
the regime PypeIt expects, where it detects a subset of the catalogue.  The band
is 7900-9100 A so that ONE list serves both cenwaves (8004-8996 A together) and
the search cannot run out of catalogue before it runs out of detector.  The arc
FWHM here is 2.1 px = 0.81 A, so 4.27 A spacing is well clear of the density that
invites misidentification (RH4 rejected 1.9 A against a 2 A FWHM for that reason).

Note `lamps = ThAr` on its own does **not** work: `HolyGrail.__init__` branches
`if 'ThAr' in self._lamps and len(self._lamps) == 1` into `run_kdtree()`, which
needs a `ThAr_patterns_*.kdtree` PypeIt 2.0.1 does not ship and cannot fetch
(`dataPaths.linelist.host is None`).  That test is exact list membership, not a
substring, so the path name `.../ThArRH3` falls through to `run_brute()`, which
needs only the list.

As on RH1, RH2 and RH4, this run is EXPECTED to end non-zero: holy-grail leaves
slices unsolved and the flat field then dies on them.  WaveCalib is written
before the flat field, which is all this stage needs.

Usage:
    python seed_rh3.py [--only 2024-10-29]
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
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
LAMPS = os.path.join(ROOT, "pypeit_test", "ThArRH3")   # PypeIt appends _lines.dat

# RH3 slices span ~793 A, inside pick_seed_slits' 400-900 A default.  Stated
# rather than assumed: RH4 lost ten minutes to a window set from another grating.
SPAN = None

TARGETS = ["2023-11-12", "2024-10-01", "2024-10-29", "2024-11-06"]


def lamp_of(path):
    """'FeAr', 'ThAr' or None, from the lamp status keywords rather than the name."""
    from astropy.io import fits
    h = fits.getheader(path)
    if str(h.get("IMTYPE")) != "ARCLAMP":
        return None
    if str(h.get("LMP0STAT")) == "1":
        return "FeAr"
    if str(h.get("LMP1STAT")) == "1":
        return "ThAr"
    return None


def retype(pfile, rawdir):
    """Make the ThAr frames the arc, and drop the FeAr ones.

    Idempotent, and it refuses to act twice: a second pass over a file whose
    FeAr rows are already commented out would find no `arc` rows to move and
    leave the setup with no arc at all.
    """
    txt = open(pfile).read()
    if "arc,tilt" in txt:
        return False, "already retyped"

    out, moved, dropped = [], 0, 0
    for line in txt.split("\n"):
        fn = line.split("|")[0].strip()
        if not fn.endswith(".fits") or not os.path.exists(os.path.join(rawdir, fn)):
            out.append(line)
            continue
        lamp = lamp_of(os.path.join(rawdir, fn))
        if lamp == "FeAr":
            out.append("# FeAr, unused for RH3: 23 catalogue lines in band are "
                       "enough to fit but not for holy-grail to identify")
            out.append("# " + line)
            dropped += 1
        elif lamp == "ThAr":
            out.append(re.sub(r"\|(\s*)tilt(\s*)\|",
                              lambda m: "|" + " " * max(1, len(m.group(1)) - 5)
                                        + "arc,tilt" + m.group(2) + "|",
                              line, count=1))
            moved += 1
        else:
            out.append(line)
    if not moved:
        return False, "no ThAr frames found -- cannot calibrate this setup"
    open(pfile, "w").write("\n".join(out))
    return True, f"{moved} ThAr frames -> arc,tilt; {dropped} FeAr frames dropped"


def set_wavelengths(pfile, sigdetect=None):
    """Rewrite the wavelengths block wholesale: holy-grail + the ThArRH3 list.

    Wholesale rather than patched, because an earlier FeAr attempt may have left
    `method = holy-grail` with the default lamps in place, and a file half
    repointed would seed against settings nobody chose.
    """
    txt = open(pfile).read()
    block = ("[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
             "    [[wavelengths]]\n"
             "        # The shipped keck_kcrm_RH3.fits covers 8599-9389 A; this\n"
             "        # data is at cenwave 8600/8400.  And FeAr's 23 in-band\n"
             "        # catalogue lines let holy-grail solve 0 of 24 slices, so\n"
             "        # the ThAr frames are the arc here -- see seed_rh3.py.\n"
             "        method = holy-grail\n"
             f"        lamps = {LAMPS}\n"
             + (f"        sigdetect = {sigdetect}\n" if sigdetect else ""))
    body = re.sub(r"\[calibrations\]\n    \[\[wavelengths\]\]\n"
                  r"(        \S.*\n|        #.*\n)+", "", txt)
    body = body.replace("[rdx]\n    spectrograph = keck_kcrm\n", block, 1)
    # Collapse runs of blank lines before comparing AND before writing.  The
    # removal regex stops at the blank line that followed the block, so a
    # rewrite reinserted the block above a blank line that was already there and
    # the file grew one blank line per call -- it therefore never compared equal
    # to itself, the "already seeded" guard never fired, and a rerun deleted a
    # good WaveCalib to recompute it.  Normalising makes this a fixed point.
    squash = lambda t: re.sub(r"\n{3,}", "\n\n", t)
    body, ref = squash(body), squash(txt)
    changed = body != ref
    if changed or body != txt:
        open(pfile, "w").write(body)
    return changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--force", action="store_true",
                    help="re-seed a night that already has a WaveCalib")
    ap.add_argument("--sigdetect", type=float, default=None,
                    help="peak-detection threshold in sigma (PypeIt default 5).  "
                         "Lower it for a night whose arcs are underexposed: "
                         "2024-10-01 took 2.3 s ThAr frames against 20 s on the "
                         "other nights, so its arc is 8x fainter (19,367 vs "
                         "151,905 ADU peak) and shows 36 lines at 5 sigma against "
                         "78.  At 3 sigma it recovers 56, comparable to what the "
                         "good nights show at 10 sigma")
    args = ap.parse_args()

    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"{LAMPS}_lines.dat missing -- run:\n"
                 f"  python pypeit_test/build_rh4_linelist.py --label RH3 "
                 f"--wmin 7900 --wmax 9100 --amp-min 300 "
                 f"--out pypeit_test/ThArRH3_lines.dat")

    for night in TARGETS:
        if args.only and args.only != night:
            continue
        sdir = os.path.join(BASE, night, "pypeit_run", "keck_kcrm_B")
        rawdir = os.path.join(ROOT, "fits", "by_night", "RH3", night)
        pf = glob.glob(os.path.join(sdir, "*.pypeit"))
        if not pf:
            print(f"{night}: no .pypeit -- run pypeit_setup first", flush=True)
            continue
        pfile = pf[0]

        ok, why = retype(pfile, rawdir)
        print(f"{night}: {why}", flush=True)
        if not ok and "already" not in why:
            continue

        changed = set_wavelengths(pfile, args.sigdetect)
        if changed:
            print(f"{night}: holy-grail + ThArRH3 line list", flush=True)

        # Don't redo a night that is already seeded.  Re-running holy-grail costs
        # ~10 minutes and, more to the point, the WaveCalib it would overwrite is
        # what build_rh3_template.py reads -- so a careless rerun silently
        # replaces good seeds with another attempt's.  Only the combination of
        # "settings already correct" and "WaveCalib present" is safe to skip on.
        done = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
        if done and not changed and not args.force:
            print(f"{night}: already seeded, not rerun (--force to redo)",
                  flush=True)
            cmd = [os.path.join(ENV, "python"), SEEDS]
            if SPAN:
                cmd += ["--span", SPAN[0], SPAN[1]]
            out = subprocess.run(cmd + [done[0]], capture_output=True, text=True)
            print(out.stdout + out.stderr, flush=True)
            continue
        # Calibrations are keyed by frame, not by parameters, so a WaveCalib built
        # under the old lamps would be reused and this would measure nothing.
        # Tilts go too: they are fit against the wavelength solution.  Arc goes
        # because the arc FRAME itself has changed -- FeAr to ThAr.
        for pat in ("WaveCalib_*.fits", "Tilts_*.fits", "Arc_*.fits",
                    "Tiltimg_*.fits"):
            for f in glob.glob(os.path.join(sdir, "Calibrations", pat)):
                os.remove(f)

        print(f"{night}: holy-grail on ThAr ...", flush=True)
        with open(os.path.join(sdir, "holygrail_thar.log"), "w") as f:
            rc = subprocess.run([os.path.join(ENV, "run_pypeit"),
                                 os.path.basename(pfile), "-c"],
                                stdout=f, stderr=subprocess.STDOUT,
                                cwd=sdir).returncode
        wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
        print(f"{night}: rc={rc} (non-zero expected), "
              f"WaveCalib {'written' if wc else 'MISSING'}", flush=True)
        if wc:
            cmd = [os.path.join(ENV, "python"), SEEDS]
            if SPAN:
                cmd += ["--span", SPAN[0], SPAN[1]]
            out = subprocess.run(cmd + [wc[0]], capture_output=True, text=True)
            print(out.stdout + out.stderr, flush=True)

    print("SEED PHASE DONE", flush=True)


if __name__ == "__main__":
    main()
