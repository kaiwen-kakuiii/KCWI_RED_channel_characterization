#!/usr/bin/env python
"""Decide whether an RM1 setup's template is fit to reduce with, before spending nights on it.

Runs calibrations for one setup with `method = full_template` against whichever
template that setup's CENTRAL WAVELENGTH is registered to, and reports how many
of the 24 slices got a usable wavelength solution.  That count is the verdict:
the unsolved slices are what crash the flat field, so a template that leaves any
of them is not usable.

## RM1 is the run where the shipped template half-works

PypeIt ships `keck_kcrm_RM1.fits` (6293.1-7733.7 A, centre 7013) and lists RM1 as
a supported setup.  An RM1 exposure is 1446 A wide, so how much of a detector has
template under it depends entirely on where the grating was set:

    cenwave 7390   74%      cenwave 6480   63%
    cenwave 6630   73%      cenwave 6300   50%
    cenwave 7510   65%      cenwave 6200   44%
                            cenwave 6130   39%

and the measured outcome is worse than that number predicts.  Gated on six
configurations, **exactly one** comes out clean.  The sharp column is not how
many slices "solved" -- it is whether the 24 slices of one setup, which are the
same grating at the same angle, agree on their DISPERSION:

    config          cenw  overlap   within 2% of median disp   median disp
    2023-12-10 C    6630    73%           24 / 24                0.7071    CLEAN
    2023-12-11 B    7510    65%           23 / 24                0.6887
    2023-12-10 B    7510    65%           20 / 24                0.6886
    2023-12-09 B    6480    63%           16 / 24                0.7099
    2024-05-02 B    6200    44%           11 / 24                0.7168
    2024-03-15 B    6300    50%            4 / 24                0.6714

**No config aborted, and none of the broken ones looked broken.**  PypeIt logged
no `Not enough useful IDs`, returned 0, and wrote all ten calibration products
INCLUDING THE FLAT FIELD on every one of them.  RH1's accidental protection --
unsolved slices crashing the flat -- does not fire here, so nothing downstream
stops.  The bad slices are visible only in their numbers: dispersions of
0.40-1.17 A/px against the ~0.71 their siblings agree on.

**A supported grating is not a supported configuration**, and 65% coverage is not
enough either.  Only 6630, at 73%, survives PypeIt's own template.

## Two things RM1 needs that earlier gratings did not

**`--span` must be widened.**  `pick_seed_slits` defaults to a 400-900 A window
for one slice, set from RH1 (~630 A) and RH2 (~764 A).  An RM1 slice spans
~1453 A.  With the default every slice is rejected as nonsense and a perfect
setup reports 0/24 -- a unit mismatch wearing the costume of bad data.  Default
here is 1200-1700.

**`shift_tol` must be raised.**  `health()` calls a slice SHIFTED when its blue
end sits more than 50 A from the median, a threshold set from RH2 (31.8-38.4 A)
and RH4 (34.0 A).  RM1's genuine slice-to-slice spread is several times that.
Measured on 2023-12-10 C -- gated 24/24 at rms 0.200 px, 70-80 lines fitted per
slit, dispersions agreeing to 0.15%, so by every other measure a clean setup:

    blue ends       5852.5 - 5991.9 A      union 139.4 A
    median          5907.6 A               max |blue - median|  84.3 A
    slices beyond the 50 A default              5
    slices beyond 100 A                         0

At the 50 A default that clean setup reports 5 SHIFTED slices.  The spread is
not misidentification: it alternates cleanly odd slice / even slice (means
5874.4 and 5949.0 A, a 74.6 A comb), which is slicer geometry.  RM1's default
here is 100.0, leaving ~16 A of headroom while still catching the 200-340 A
misidentifications the tolerance exists for.

## The lamps are PypeIt's own, and that is measured rather than assumed

RH3 and RH4 had to build ThAr catalogues because the default FeI/ArI/ArII list
was too thin in band.  RM1 does not: on 2023-12-10 C the default catalogue fitted
**70-80 lines per slit** at rms 0.200 px.  RH3's rule is to count catalogue lines
against the method that will consume them, and full_template -- which is handed
the identification and only fits -- is very well served here.

Usage:
    python gate_rm1_template.py --night 2023-12-11 --setup keck_kcrm_C
    python gate_rm1_template.py --night 2024-03-15 --setup keck_kcrm_B --template <path>
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
RUNS = os.path.join(ROOT, "RM1 pypeit run")

SHIPPED = "keck_kcrm_RM1.fits"      # resolved through PypeIt's own data cache

# Which template serves which central wavelength.  SHIPPED where it gates clean;
# one of ours where it does not.  Keyed by cenwave rounded to 10 A, because RM1's
# settings sit 100-130 A apart and a coarser round would merge them.
#
# The blue group shares ONE template rather than taking one each: an RM1 exposure
# is 1446 A wide, so a template seeded at 6300 covers 5577-7023 A and sits under
# 88-100% of the 6130/6200/6300/6480 detectors -- far above the ~50% where the
# shipped one fails.  Seven templates for seven settings would be busywork.
REGISTRY = {
    6130: os.path.join(HERE, "keck_kcrm_RM1_6130_small.fits"),
    6200: os.path.join(HERE, "keck_kcrm_RM1_6300.fits"),
    6300: os.path.join(HERE, "keck_kcrm_RM1_6300.fits"),
    6480: os.path.join(HERE, "keck_kcrm_RM1_6300.fits"),
    6630: SHIPPED,
    7390: os.path.join(HERE, "keck_kcrm_RM1_7390_small.fits"),
    7510: os.path.join(HERE, "keck_kcrm_RM1_7510.fits"),
}

# Which setups sit at which central wavelength, so the gate can say whether it is
# testing generalisation or reproducing the arc it was seeded from.
SETUPS = {
    6130: ["2024-04-01/keck_kcrm_A"],
    6200: ["2024-05-02/keck_kcrm_B"],
    6300: ["2024-03-15/keck_kcrm_B"],
    6480: ["2023-12-09/keck_kcrm_B"],
    6630: ["2023-12-10/keck_kcrm_C", "2023-12-11/keck_kcrm_C"],
    7390: ["2024-04-01/keck_kcrm_B"],
    7510: ["2023-12-10/keck_kcrm_B", "2023-12-11/keck_kcrm_B"],
}

SPAN_DEFAULT = (1200.0, 1700.0)
# SHIFTED tolerance is a property of the grating AND THE SLICER, not just the
# grating.  RH3 section 5.2 made it a parameter after measuring 44.6 A on RH3
# against 31.8-38.4 A on RH2; RM1 shows the slicer moves it again.  Measured on
# setups whose dispersions agree to 0.3% and whose worst slice fits 97 lines at
# rms 0.118 -- i.e. correct solutions, not misidentifications:
#
#     6630 Large 2x2   blue-end union 139 A   max |blue - median|   84.3 A
#     7510 Large 2x2                  129 A                         78.6 A
#     6130 Small 1x1                  170 A                        105.8 A
#
# A misidentification moves a solution by about a line spacing and shows up in
# the DISPERSION.  6130's most deviant slice reads -0.19% against its siblings,
# which is exactly what 7510's most deviant slice reads.  Both are geometry.
SHIFT_TOL_BY_DECKER = {"Large": 100.0, "Small": 120.0}
SHIFT_TOL_DEFAULT = 100.0


def shift_tol_for(decker):
    """Tolerance for this slicer, with the grating's 200-340 A failures still caught."""
    return SHIFT_TOL_BY_DECKER.get(decker, SHIFT_TOL_DEFAULT)
HALF_SPAN = 728.0


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
                    help="physically possible span of ONE RM1 slice, in A "
                         "(measured ~1453).  The pick_seed_slits default of "
                         "400-900 rejects every RM1 slice as nonsense")
    ap.add_argument("--shift-tol", type=float, default=None,
                    help="A a slice's blue end may sit from the median before "
                         "it counts as SHIFTED.  RM1's genuine spread is 139 A "
                         "across the 24 slices (84 A from the median), so the "
                         "50 A default flags five slices of a setup that is "
                         "otherwise perfect at 24/24 and rms 0.200 px")
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
        sys.exit(f"template missing: {tmpl_name} -- run build_rm1_template.py "
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
    if cov < 0.55:
        print(f"   *** below the ~50-55% where the shipped template was measured "
              f"to fail on RM1 and RH3 -- expect this to fail ***")
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
    print(f"   lamps    PypeIt's default FeI/ArI/ArII.  Measured on "
          f"2023-12-10 C: 70-80 lines fitted per slit, so unlike RH3/RH4 no "
          f"custom catalogue is needed")

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
    good = rank(wc[0], quiet=True, span=span)

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
              f"That is what build_rm1_template.py bootstraps from (RH3 section 4).")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
