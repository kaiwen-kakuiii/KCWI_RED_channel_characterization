#!/usr/bin/env python
"""Bootstrap wavelength seeds for RH4, which PypeIt can neither template nor identify.

RH1 and RH2 were missing a template.  **RH4 is missing the line list as well**,
and that is the harder half.

    catalogue lines in RH4's 9419-10402  FeI  0    ArI  4    ArII  2   (6 total)
    catalogue lines in RH2's 6400-7600   FeI 25    ArI 40    ArII 13   (78)

`keck_kcwi.py` assigns `lamps = ['FeI', 'ArI', 'ArII']`, `FeI_lines.dat` stops at
9002.0 A, and RH4 is centred at 9950.  A template would not help: `full_template`
uses it to decide which *catalogue* line each feature is and then fits to
catalogue wavelengths, so with six lines in band -- four of them in the bluest
370 A -- there is nothing to fit.

The lamp is fine -- measured on the raw frames, FeAr puts 114 lines over 100
sigma on the RH4 detector, as many as on RH2.  PypeIt's curated lists simply stop
before 9500 A.

**Fix: calibrate off the ThAr frames.**  KCWI takes both lamps and types FeAr
`arc`, ThAr `tilt` (`check_frame_type`: "PypeIt is only setup to wavelength
calibrate using the FeAr lamp").  `ThAr_lines.dat` carries 854 lines in 9200-10800 A,
811 of them NIST-flagged.  So this script retypes the three ThAr frames `arc,tilt`, comments
out the three FeAr frames, and points `lamps` at the band-limited list from
`build_rh4_linelist.py`.

Note `lamps = ThAr` on its own does **not** work: `HolyGrail.__init__` branches
`if 'ThAr' in self._lamps and len(self._lamps) == 1` into `run_kdtree()`, which
needs a `ThAr_patterns_*.kdtree` file PypeIt 2.0.1 does not ship and cannot fetch
here (`dataPaths.linelist.host is None`).  An absolute path as the lamp name
falls through to `run_brute()` and is accepted, because `lamps` is not validated
against a fixed list and `get_file_path` returns any path that already resolves.

As on RH1 and RH2, this run is EXPECTED to end with a non-zero rc: holy-grail
leaves slices unsolved and the flat field then dies on them
(`flatfield.py illum_profile_spectral` takes np.min of an empty array).  WaveCalib
is written before the flat field, which is all this stage needs.
"""
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
SEEDS = os.path.join(ROOT, "pypeit_test", "pick_seed_slits.py")
LAMPS = os.path.join(ROOT, "pypeit_test", "ThArRH4")

TARGETS = [("2024-12-28", "keck_kcrm_B")]


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
            out.append("# FeAr, unused for RH4: PypeIt's FeI/ArI/ArII lists hold "
                       "4 lines above 9550 A")
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


for night, setup in TARGETS:
    sdir = os.path.join(BASE, night, "pypeit_run", setup)
    rawdir = os.path.join(RAW, night)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        print(f"{night}/{setup}: no .pypeit -- run pypeit_setup first", flush=True)
        continue
    pfile = pfile[0]

    if not os.path.exists(LAMPS + "_lines.dat"):
        sys.exit(f"{LAMPS}_lines.dat missing -- run build_rh4_linelist.py first")

    ok, why = retype(pfile, rawdir)
    print(f"{night}/{setup}: {why}", flush=True)
    if not ok and "already" not in why:
        continue

    txt = open(pfile).read()
    if "lamps" not in txt:
        txt = txt.replace(
            "[rdx]\n    spectrograph = keck_kcrm\n",
            "[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n"
            "    [[wavelengths]]\n        method = holy-grail\n"
            f"        lamps = {LAMPS}\n", 1)
        open(pfile, "w").write(txt)
        print(f"{night}/{setup}: set holy-grail + ThArRH4 line list", flush=True)
        # Calibrations are keyed by frame, not by parameters, so anything built
        # under the old lamp list would be reused rather than rebuilt.
        if os.path.isdir(os.path.join(sdir, "Calibrations")):
            shutil.rmtree(os.path.join(sdir, "Calibrations"))
            print(f"{night}/{setup}: cleared pre-ThAr Calibrations", flush=True)

    wc = glob.glob(os.path.join(sdir, "Calibrations", "WaveCalib_*.fits"))
    if wc:
        print(f"{night}/{setup}: WaveCalib already present, not rerun", flush=True)
    else:
        print(f"{night}/{setup}: holy-grail on ThAr ...", flush=True)
        log = os.path.join(sdir, "holygrail_thar.log")
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
