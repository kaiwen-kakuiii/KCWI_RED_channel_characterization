#!/usr/bin/env python
"""Decide whether an RM2 setup's template is fit to reduce with, before spending nights on it.

Same job and same verdict as gate_rm1_template.py: run calibrations for one
setup with `method = full_template` against whichever template that setup's
CENTRAL WAVELENGTH is registered to, then judge the 24 slices.  The sharp column
is not how many slices "solved" -- RM1 section 3 -- it is whether the 24 slices
of one setup, which are the same grating at the same angle, agree on their
DISPERSION.

## What is different about RM2

**PypeIt ships `keck_kcrm_RM2.fits` and lists RM2 as supported.**  It covers
8608-10595 A, centre 9601, binspec 2, 0.9659 A/px.  An RM2 exposure is
2064 x 0.9659 = 1994 A wide, so coverage per central wavelength in this run:

    cenwave 9850   Medium 2x2    87%
    cenwave 8950   Medium 2x2    67%
    cenwave 8900   Medium 2x2    65%
    cenwave 8850   Small  1x1    62%
    cenwave 7750   Small  1x1     7%

RM1 measured 65% NOT enough and 74% not enough either when the slicer differed
from the template's, against 73% passing on Medium/2x2.  So expect 9850 to pass,
8950/8900 to be genuinely marginal, and the two Small configs to need templates
bootstrapped from their own good slices (RH3 section 4, build_rm1_template.py).

**The slicers are Medium and Small.  There is no Large config in RM2**, so
`SHIFT_TOL_BY_DECKER`'s RM1 values ("Large" 100, "Small" 120) do not transfer by
name.  The tolerances here start PROVISIONAL and are re-measured off the first
setup that gates clean -- the rule RM1 section 5.2 states: read the number that
separates geometry from error off setups already known good, never import it.

**A slice spans ~1994 A, not RM1's ~1453.**  `pick_seed_slits` defaults to a
400-900 A window and would reject every RM2 slice as nonsense; RM1 already had
to widen it to 1200-1700, and that window is just as wrong here.  Default below
is 1700-2300.

**The arc is FeAr and the tilt frames are ThAr**, checked on all nine setups --
the same assignment RM1 had, so PypeIt's default FeI/ArI/ArII catalogue is what
full_template is handed.  Verify the fitted line count on the first clean gate:
RM2 sits 1000-3000 A redder than RM1, where that catalogue thins out, and RH3/RH4
had to build their own at 8400-8600 A.

Usage:
    python gate_rm2_template.py --night 2023-09-23 --setup keck_kcrm_B
    python gate_rm2_template.py --night 2024-05-09 --setup keck_kcrm_A --template <path>
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank, health              # noqa: E402

ENV = "/opt/miniconda3/envs/pypeit/bin"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RM2 pypeit run")

SHIPPED = "keck_kcrm_RM2.fits"      # resolved through PypeIt's own data cache

# Which template serves which central wavelength.  SHIPPED where it gates clean;
# one of ours where it does not.  Keyed by cenwave rounded to 10 A.
#
# Every entry starts at SHIPPED because RM2 is a supported setup and the shipped
# template has to be given its chance at each angle -- that is what the gate is
# for.  Entries move to a bootstrapped template only after a gate has FAILED at
# that cenwave and build_rm2_template.py has seeded one; the registry is then the
# single source of truth the driver checks the .pypeit against.
REGISTRY = {
    7750: os.path.join(HERE, "keck_kcrm_RM2_7750_small.fits"),
    8850: os.path.join(HERE, "keck_kcrm_RM2_8900.fits"),
    8900: os.path.join(HERE, "keck_kcrm_RM2_8900.fits"),
    8950: os.path.join(HERE, "keck_kcrm_RM2_8900.fits"),
    9850: SHIPPED,
}

# Measured with the shipped template, before any of ours existed.  This is the
# table the registry above is derived from, kept here so a reader can see why
# each entry is what it is:
#
#   config          cenw  slicer  overlap  solved  what happened
#   2023-09-23 B    9850  Medium    87%    24/24   CLEAN, dispersions agree 0.4%
#   2025-01-01 B    8950  Medium    67%     5/24   flat crashed
#   2024-06-11 B    8950  Medium    67%     2/24   flat crashed
#   2024-01-04 B    8900  Medium    65%    12/24   flat crashed; 8 of the 12 carry
#                                                  disp 0.66-0.83 vs the true 0.975
#   2024-06-10 A    8850  Small     62%     0/24   nothing solved at all
#
# So the Medium threshold sits between 67% and 87%, and 8900/8950 move to
# keck_kcrm_RM2_8900.fits (7871-10077 A, seeded from 2024-01-04's two good
# slices) which covers 100% of both detectors.  9850 keeps PypeIt's own.
#
# 8850 is Small/1x1 and takes the SAME Medium-seeded template, which RH2.md says
# should not work -- "a template cannot serve two slicers".  Measured on
# 2024-06-10 A, it does: 24/24 solved, rms median 0.259 px, 19-23 lines per
# slice, and the 24 dispersions agree to 0.24% around 0.4885 A/px, which is half
# the Medium 0.9752 to within 0.2%.  The RH2 rule was measured on a
# high-dispersion RH grating, where a resolution change blends neighbouring
# lines; RM2's FeAr lines out here sit 5-9 A apart against 0.49 A/px sampling, so
# there is nothing to blend.  The rule's real variable is line spacing versus the
# resolution difference, not the slicer name.
#
# 7750 is the one setting nothing reached: 7% of its detector under the shipped
# template, 44% under the Medium one, and gated against the latter its 24
# dispersions scatter 0.4829-0.5348 with only 9 of 19 inside 2% -- a partly wrong
# solution, not a clean failure.  It is served by keck_kcrm_RM2_7750_small.fits,
# seeded by holy-grail (holy_seed_rm2.py) from 2024-05-09 itself, which solved
# 24/24 at 39-52 lines per slice and rms median 0.111 px because at 7750 the FeAr
# catalogue is still rich.  Line density is the whole story of this run.

# Which setups sit at which central wavelength, so the gate can say whether it is
# testing generalisation or reproducing the arc it was seeded from.  Note 8850:
# two nights at one cenwave but DIFFERENT STARS (feige110 / feige 34), and 8950
# has three nights -- the strongest reproducibility test in the project so far.
SETUPS = {
    7750: ["2024-05-09/keck_kcrm_A"],
    8850: ["2024-06-10/keck_kcrm_A", "2024-12-24/keck_kcrm_A"],
    8900: ["2024-01-04/keck_kcrm_B", "2024-04-30/keck_kcrm_B"],
    8950: ["2024-06-11/keck_kcrm_B", "2025-01-01/keck_kcrm_B",
           "2025-01-02/keck_kcrm_B"],
    9850: ["2023-09-23/keck_kcrm_B"],
}

SPAN_DEFAULT = (1700.0, 2300.0)
# SHIFTED tolerance is a property of the grating AND THE SLICER, not just the
# grating (RH3 section 5.2, RM1 section 5.2).  RM2 has no Large config, so RM1's
# values cannot be carried over by decker name.
#
# Both are now MEASURED on setups gated 24/24 with their 24 dispersions agreeing
# to better than 0.5%, i.e. on solutions known to be right:
#
#     2023-09-23 B   9850   Medium 2x2   blue-end union 206.7 A   max |blue-med| 128.4 A
#     2024-06-10 A   8850   Small  1x1   blue-end union 247.0 A   max |blue-med| 150.9 A
#
# so 150 and 180 stand with 22 and 29 A of headroom.  The Small value being the
# larger repeats RM1's finding that the slicer moves this number, not only the
# grating.  A misidentification at these dispersions moves a solution by several
# hundred A, which both values still catch.
SHIFT_TOL_BY_DECKER = {"Medium": 150.0, "Small": 180.0}

# Lines a slice must fit to count toward the "strict seed criteria" line of the
# report.  pick_seed_slits' 20 comes from gratings where FeAr fits 70-80 lines
# per slice; RM2 fits 10-14 at 9850 and 22-28 at 8900/8950, so 20 would report
# zero seeds on the one setup that gated 24/24.  See pick_seed_slits.py.
NLINE_MIN = 8
SHIFT_TOL_DEFAULT = 150.0


def shift_tol_for(decker):
    """Tolerance for this slicer, with the grating's several-hundred-A failures still caught."""
    return SHIFT_TOL_BY_DECKER.get(decker, SHIFT_TOL_DEFAULT)


HALF_SPAN = 997.0   # half an RM2 exposure: 2064 px x 0.9659 A/px / 2


def resolve(tmpl):
    """Absolute path to a template, for READING it here.

    Not for writing into a .pypeit -- see `reid_value`.  PypeIt caches a
    downloaded arxiv template as an extensionless blob
    (`~/.cache/pypeit/download/url/<hash>/contents`), which is fine to open
    directly but useless as a parameter value.
    """
    if os.path.isabs(tmpl) or os.path.exists(tmpl):
        return os.path.abspath(tmpl)
    from pypeit import dataPaths
    return str(dataPaths.reid_arxiv.get_file_path(tmpl))


def reid_value(tmpl):
    """What to write as `reid_arxiv` in a .pypeit.

    Ours go in by absolute path.  PypeIt's own goes in BY NAME, so PypeIt
    resolves it through its own data machinery.  Writing the resolved cache path
    instead kills the run:

        reid_arxiv = ~/.cache/pypeit/download/url/412d0c.../contents
        IORegistryError: No reader defined for format '' and class 'Table'

    -- astropy infers the format from the extension, and the cache blob has
    none.  This cost two 6630 gates that had already passed at 24/24 with
    exactly the same template.
    """
    return tmpl if tmpl == SHIPPED else os.path.abspath(tmpl)


def seeded_from(tmpl):
    """Setups whose slices seeded this template, read from its header.

    Read from the FILE, never inferred from which WaveCalibs currently look good.
    Once a setup has been gated its WaveCalib is a full_template solution and
    passes rank() too, so every setup looks like a seed source and a reproduction
    gets reported as independent -- the error RH3 section 9 documents.
    """
    try:
        h = fits.getheader(tmpl)
    except Exception:
        return None
    if "SEEDNITE" not in h:
        return None
    return [s for s in str(h["SEEDNITE"]).split(",") if s]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", required=True)
    ap.add_argument("--setup", default="keck_kcrm_B")
    ap.add_argument("--template", default=None,
                    help="override the registry, e.g. to test the shipped "
                         "template against a setup it is not registered for")
    ap.add_argument("--min-frac", type=float, default=1.0,
                    help="fraction of slices that must solve to pass")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=list(SPAN_DEFAULT),
                    help="physically possible span of ONE RM2 slice, in A "
                         "(~1994 = 2064 px x 0.9659).  The pick_seed_slits "
                         "default of 400-900 rejects every RM2 slice as nonsense")
    ap.add_argument("--shift-tol", type=float, default=None,
                    help="A a slice's blue end may sit from the median before "
                         "it counts as SHIFTED.  RM2's own spread has not been "
                         "measured yet -- the defaults are scaled from RM1 by "
                         "dispersion and must be replaced by a measurement off "
                         "the first clean gate")
    ap.add_argument("--keep-calibs", action="store_true",
                    help="do not delete the existing Calibrations/ first.  Off "
                         "by default: PypeIt keys calibrations by frame, not by "
                         "parameters, so the previous products would be reused "
                         "and the test would measure nothing")
    args = ap.parse_args()

    sdir = os.path.join(RUNS, args.night, "pypeit_run", args.setup)
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        sys.exit(f"no .pypeit in {sdir}")
    pfile = pfile[0]
    txt = open(pfile).read()

    mc = re.search(r"cenwave:\s*([0-9.]+)", txt)
    cw = int(round(float(mc.group(1)) / 10.0) * 10) if mc else None
    md = re.search(r"decker:\s*(\S+)", txt)
    decker = md.group(1) if md else "?"

    tmpl_name = args.template or REGISTRY.get(cw)
    if tmpl_name is None:
        sys.exit(f"cenwave {cw} is not in REGISTRY -- add it, or pass --template")
    if tmpl_name != SHIPPED and not os.path.exists(tmpl_name):
        sys.exit(f"template missing: {tmpl_name} -- run build_rm2_template.py "
                 f"--cenwave {cw}" + ("" if decker.lower() == "large"
                                      else f" --decker {decker}") + " first")
    tmpl = resolve(tmpl_name)

    tw = np.asarray(fits.getdata(tmpl)["wave"], float).ravel()
    tw = tw[tw > 0]
    lo, hi = cw - HALF_SPAN, cw + HALF_SPAN
    cov = max(0.0, min(hi, tw.max()) - max(lo, tw.min())) / (hi - lo)

    print(f"{args.night}/{args.setup}: decker {decker}, cenwave {cw}")
    print(f"   template {os.path.basename(tmpl_name)}"
          + ("   (PypeIt's own)" if tmpl_name == SHIPPED else "   (built here)"))
    print(f"   setup covers ~{lo:.0f}-{hi:.0f} A, template covers "
          f"{tw.min():.0f}-{tw.max():.0f} A -> {100*cov:.0f}% of the detector "
          f"has template under it")
    if cov < 0.70:
        print(f"   *** RM1 measured 65% NOT enough and 74% not enough across a "
              f"slicer change; below ~70% expect this to fail ***")
    src = seeded_from(tmpl)
    if src is None:
        print("   provenance: cannot be determined from the file "
              + ("(PypeIt ships no SEEDNITE header)" if tmpl_name == SHIPPED
                 else "-- rebuild it to stamp SEEDNITE"))
    elif f"{args.night}/{args.setup}" in src:
        print(f"   provenance: this setup SEEDED the template ({','.join(src)}) "
              f"-- this gate is REPRODUCTION, not generalisation")
    else:
        print(f"   provenance: seeded from {','.join(src)} -- this gate is an "
              f"INDEPENDENT test of whether the template generalises")
    print(f"   lamps    arc = FeAr, tilt = ThAr (checked on all nine RM2 "
          f"setups); catalogue is PypeIt's default FeI/ArI/ArII.  RM2 sits "
          f"1000-3000 A redder than RM1, where that list thins -- read the "
          f"fitted line count below, do not assume it")

    # Rebuild the whole parameter header rather than patching pieces of it, so a
    # `lamps` or `method` line left over from an earlier attempt cannot survive
    # into the gated settings.  Everything from the setup block onward is data
    # and is never touched.
    #
    # `exclude_regions` is READ BACK OUT and re-emitted.  It is a hand-made
    # decision about slit tracing (2024-04-01 B needs one -- see RM1.md), and an
    # earlier version of this rewrite silently dropped it, which put the spurious
    # 25th slice straight back.  A gate that discards the fix it is testing is
    # worse than no gate.
    mset = re.search(r"\n(# Setup\n)?setup read\n", txt)
    if not mset:
        sys.exit(f"{pfile}: no 'setup read' block found")
    mex = re.search(r"^\s*exclude_regions\s*=\s*(\S+)\s*$",
                    txt[:mset.start()], re.M)

    head = ("# Auto-generated PypeIt input file using PypeIt version: 2.0.1\n"
            "\n# User-defined execution parameters\n"
            "[rdx]\n    spectrograph = keck_kcrm\n\n[calibrations]\n")
    if mex:
        head += f"    [[slitedges]]\n        exclude_regions = {mex.group(1)}\n"
    head += ("    [[wavelengths]]\n        method = full_template\n"
             f"        reid_arxiv = {reid_value(tmpl_name)}\n")
    open(pfile, "w").write(head + txt[mset.start():])
    print(f"patched {os.path.basename(pfile)}"
          + (f"   (kept exclude_regions = {mex.group(1)})" if mex else ""))

    calib = os.path.join(sdir, "Calibrations")
    if os.path.isdir(calib) and not args.keep_calibs:
        shutil.rmtree(calib)
        print("removed the previous Calibrations/ so it cannot be reused")

    log = os.path.join(sdir, "gate.log")
    print(f"running calibrations ... (log: {log})", flush=True)
    with open(log, "w") as f:
        rc = subprocess.run([os.path.join(ENV, "run_pypeit"),
                             os.path.basename(pfile), "-c"],
                            stdout=f, stderr=subprocess.STDOUT, cwd=sdir).returncode

    wc = glob.glob(os.path.join(calib, "WaveCalib_*.fits"))
    if not wc:
        sys.exit(f"FAIL: run_pypeit rc={rc} and no WaveCalib was written -- see {log}")

    # The wavelength solution is not the whole calibration, and judging the gate
    # on WaveCalib alone reported PASS on three configs whose FLAT FIELD was
    # never written.  Cause, that time: `run_pypeit` died writing a QA png into a
    # directory that had been deleted underneath it, part-way through building
    # the flat -- rc was 1, WaveCalib was already on disk, and the gate said
    # "PASS 24/24".  That is precisely the failure mode this gate exists to
    # catch, so it now checks the products and the exit code as well.
    missing = [k for k in ("Arc", "Edges", "Slits", "Flat", "Tilts", "WaveCalib")
               if not glob.glob(os.path.join(calib, f"{k}_*"))]
    if missing or rc != 0:
        print(f"\n*** run_pypeit returned {rc}"
              + (f" and did not write: {', '.join(missing)}" if missing else "")
              + f" -- see {log}")
    if missing:
        sys.exit(f"FAIL: the calibration is incomplete regardless of how the "
                 f"wavelength solution looks.  Missing {', '.join(missing)}.")

    span = (float(args.span[0]), float(args.span[1]))
    stol = args.shift_tol if args.shift_tol is not None else shift_tol_for(decker)
    n, solved, flagged, rms, shifted = health(wc[0], span=span, shift_tol=stol)
    ok = [s for s in solved if s not in flagged and s not in shifted]
    good = rank(wc[0], quiet=True, span=span, nline_min=NLINE_MIN)

    print(f"\nrun_pypeit -c returned {rc}")
    print(f"slices solved, unflagged and unshifted: {len(ok)} / {n}")
    if flagged:
        print(f"   {len(flagged)} flagged by PypeIt: {flagged}")
    if shifted:
        print(f"   {len(shifted)} blue end >{stol:.0f} A from median: {shifted}")
    if len(rms):
        print(f"   rms median {np.median(rms):.3f} px, worst {np.max(rms):.3f} px")
    print(f"   {len(good)} / {n} also meet the strict seed criteria")

    verdict = ("PASS" if (len(ok) >= args.min_frac * n and rc == 0)
               else "FAIL")
    print(f"\n{verdict}: {len(ok)}/{n} usable slices "
          f"(needed {args.min_frac*n:.0f}), run_pypeit rc={rc}, "
          f"all six calibration products present")
    if verdict == "FAIL" and len(good):
        print(f"   NOT necessarily a dead end: {len(good)} slices still meet the "
              f"strict seed criteria, spanning "
              f"{min(g['wmin'] for g in good):.0f}-{max(g['wmax'] for g in good):.0f} A.  "
              f"That is what build_rm2_template.py bootstraps from (RH3 section 4).")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
