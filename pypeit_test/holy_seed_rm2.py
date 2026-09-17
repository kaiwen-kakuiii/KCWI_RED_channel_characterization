#!/usr/bin/env python
"""Run holy-grail on one RM2 setup to get seed slices where no template reaches.

The bootstrap of RH1_PROCEDURE.md, made a script because RM2 needs it twice.

RM2's two Small/1x1 configurations have nothing to start from.  PypeIt's shipped
`keck_kcrm_RM2.fits` covers 8608-10595 A, which is 62% of an 8850 detector and 7%
of a 7750 one, and it is a Medium-resolution template besides; the Medium
template built here (`keck_kcrm_RM2_8900.fits`, 7871-10077 A) covers the 8850
detector geometrically but was seeded at the wrong slicer, and RH2.md's finding
is that a template holds as single features the blends another slicer resolves.
So the seeds have to come from the data itself, with no template at all.

holy-grail identifies lines from scratch.  It is expected to solve FEW slices --
1 of 24 on RH1 -- and to leave the rest silently wrong, which is exactly why its
output goes through `pick_seed_slits.rank()` and not into a reduction:

    rms      < 0.5 binned px
    nlines  >= 8              (RM2's FeAr density; 20 is an RM1/RH number)
    span     = 1700-2300 A    (one RM2 slice is ~1994 A)
    monotonic dispersion, and d(lambda)/dpix within 25% of the median

A single surviving slice is enough to seed a template.  What must NOT happen is
reducing science against the holy-grail solution itself -- on RM1 a wrongly
solved setup wrote all ten calibration products with rc=0 and nothing downstream
stopped.

The flat field is expected to crash here (`ValueError: zero-size array`, from
unsolved slices), and that is FINE: this script only needs WaveCalib.  It reports
the seeds and stops; nothing is registered and no science is touched.

Usage:
    python holy_seed_rm2.py --night 2024-06-10 --setup keck_kcrm_A
    python holy_seed_rm2.py --night 2024-05-09 --setup keck_kcrm_A --sigdetect 3
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank, health              # noqa: E402

ENV = "/opt/miniconda3/envs/pypeit/bin"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RM2 pypeit run")

SPAN_DEFAULT = (1700.0, 2300.0)
NLINE_MIN = 8
SHIFT_TOL = 180.0          # Small/1x1; see gate_rm2_template.py


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", required=True)
    ap.add_argument("--setup", default="keck_kcrm_A")
    ap.add_argument("--sigdetect", type=float, default=None,
                    help="lower it (5 -> 3) to hand holy-grail more lines.  RM2 "
                         "sits where the FeAr catalogue thins, and holy-grail "
                         "needs enough lines to IDENTIFY, which is a stricter "
                         "requirement than enough to fit")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=list(SPAN_DEFAULT))
    ap.add_argument("--nline-min", type=int, default=NLINE_MIN)
    ap.add_argument("--keep-calibs", action="store_true")
    args = ap.parse_args()

    sdir = os.path.join(RUNS, args.night, "pypeit_run", args.setup)
    pf = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pf:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pf[0]
    txt = open(pfile).read()

    md = re.search(r"decker:\s*(\S+)", txt)
    mc = re.search(r"cenwave:\s*([0-9.]+)", txt)
    print(f"{args.night}/{args.setup}: decker {md.group(1) if md else '?'}, "
          f"cenwave {float(mc.group(1)) if mc else 0:.0f}, method holy-grail "
          f"(no reid_arxiv)")

    # Rebuild the parameter header, dropping reid_arxiv entirely -- leaving it in
    # with method = holy_grail is a contradiction PypeIt does not complain about.
    # exclude_regions is carried through for the reason gate_rm2_template.py
    # carries it: it is a hand-made slit-tracing decision, and dropping it here
    # would seed from a different slit set than the one the run will use.
    mset = re.search(r"\n(# Setup\n)?setup read\n", txt)
    if not mset:
        sys.exit(f"{pfile}: no 'setup read' block found")
    mex = re.search(r"^\s*exclude_regions\s*=\s*(\S+)\s*$", txt[:mset.start()], re.M)

    head = ("# Auto-generated PypeIt input file using PypeIt version: 2.0.1\n"
            "\n# User-defined execution parameters\n"
            "[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n")
    if mex:
        head += f"    [[slitedges]]\n        exclude_regions = {mex.group(1)}\n"
    # `holy-grail`, with a HYPHEN.  PypeIt's own docs and this project's notes
    # write it "holy_grail" in prose, and the underscore is rejected --
    # `ValueError: Input value for method invalid: holy_grail.  Options are:
    # holy-grail, identify, reidentify, echelle, full_template` -- after the
    # calibration run has already started, so the cost is a wasted pass.
    head += "    [[wavelengths]]\n        method = holy-grail\n"
    if args.sigdetect is not None:
        head += f"        sigdetect = {args.sigdetect}\n"
    open(pfile, "w").write(head + txt[mset.start():])
    print(f"patched {os.path.basename(pfile)}"
          + (f"   (kept exclude_regions = {mex.group(1)})" if mex else ""))

    calib = os.path.join(sdir, "Calibrations")
    if os.path.isdir(calib) and not args.keep_calibs:
        shutil.rmtree(calib)
        print("removed the previous Calibrations/ -- PypeIt keys calibrations by "
              "frame, not by parameters, so it would otherwise be reused")

    log = os.path.join(sdir, "holy.log")
    print(f"running calibrations ... (log: {log})", flush=True)
    with open(log, "w") as f:
        rc = subprocess.run([os.path.join(ENV, "run_pypeit"),
                             os.path.basename(pfile), "-c"],
                            stdout=f, stderr=subprocess.STDOUT, cwd=sdir).returncode

    wc = glob.glob(os.path.join(calib, "WaveCalib_*.fits"))
    if not wc:
        sys.exit(f"no WaveCalib written (rc={rc}) -- see {log}.  holy-grail found "
                 f"nothing to work with; try --sigdetect 3")

    span = (float(args.span[0]), float(args.span[1]))
    n, solved, flagged, rms, shifted = health(wc[0], span=span,
                                              shift_tol=SHIFT_TOL)
    good = rank(wc[0], quiet=False, span=span, nline_min=args.nline_min)

    print(f"\nrun_pypeit -c returned {rc}"
          + ("   (a crashed flat is EXPECTED here and does not matter -- only "
             "WaveCalib is wanted)" if rc else ""))
    print(f"slices with a physically possible solution: {len(solved)} / {n}")
    if len(rms):
        print(f"   rms median {np.median(rms):.3f} px, worst {np.max(rms):.3f} px")
    print(f"{len(good)} / {n} meet the seed criteria "
          f"(rms<0.5, nlines>={args.nline_min}, span {span[0]:.0f}-{span[1]:.0f} A)")
    if good:
        print(f"   seeds span {min(g['wmin'] for g in good):.1f}-"
              f"{max(g['wmax'] for g in good):.1f} A, dispersions "
              f"{min(g['disp'] for g in good):.4f}-{max(g['disp'] for g in good):.4f}")
        print(f"\nNext: build_rm2_template.py --cenwave <cw> --decker "
              f"{md.group(1) if md else 'Small'} --only {args.night}")
    else:
        print("   nothing usable.  Options, in order of cost: --sigdetect 3; "
              "then the ThAr frames (PypeIt types them `tilt`) with a catalogue "
              "built the way build_rh4_linelist.py builds one.")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
