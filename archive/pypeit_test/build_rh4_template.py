#!/usr/bin/env python
"""Build the wavelength template PypeIt does not ship for the RH4 grating.

Same hole as RH1 and RH2 -- `keck_kcwi.py::config_specific_par` assigns
`reid_arxiv` for RL, RM1, RM2 and RH3 only -- but RH4 reaches it with two
differences that change the shape of the job.

**1. The template is useless without the RH4 line list.**  `full_template` uses
the template only to decide which *catalogue* line each feature is, then fits to
catalogue wavelengths.  PypeIt's KCWI catalogue (`FeI`, `ArI`, `ArII`) holds four
lines above 9550 A, so a template alone changes nothing.  RH4 calibrates against
`ThArRH4_lines.dat` (see build_rh4_linelist.py), and **every run that uses this
template must also set `lamps` to that list** -- they are one unit, not two
independent settings.

**2. One night, one configuration, so one seed and no stitch.**  RH1 needed a
chain across six central wavelengths and RH2 needed two chains, one per slicer.
RH4 is a single Large-slicer config at cenwave 9950, so the template is built
from one setup and comes out one exposure wide -- the shape PypeIt ships, with no
stitch joint and no alignment ambiguity (RH1_PROCEDURE.md problem 3).  That is
the good case; it is not a compromise.

It also means the honest gate RH2 used -- build from one night, gate on a
*different* night -- is not available here.  See gate_rh4_template.py, which says
so rather than implying a strength the data cannot support.

Usage:
    python build_rh4_template.py
    python build_rh4_template.py --max-seeds 2      # if a second config appears
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
from astropy.table import Table

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank                       # noqa: E402
from pypeit.core.wavecal import templates              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RH4 pypeit run")
OUT = os.path.join(HERE, "keck_kcrm_RH4.fits")

# RH4's nights are all Large so far.  Kept explicit rather than assumed: the
# slicer sets spectral resolution, and RH1 lost a day to a template built at one
# slicer being run against another (RH1_PROCEDURE.md problem 2).
DECKER = "Large"


def setup_meta(wavecalib_path):
    """(spectral binning, decker) of the setup a WaveCalib belongs to.

    Read from the .pypeit setup block rather than a FITS header, so it describes
    the setup PypeIt actually grouped rather than one frame of it.
    """
    sdir = os.path.dirname(os.path.dirname(wavecalib_path))
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        raise RuntimeError(f"no .pypeit beside {wavecalib_path}")
    txt = open(pfile[0]).read()
    mb = re.search(r"binning:\s*(\d+),(\d+)", txt)
    md = re.search(r"decker:\s*(\S+)", txt)
    if not (mb and md):
        raise RuntimeError(f"cannot read binning/decker for {wavecalib_path}")
    # KCWI writes binning as spatial,spectral
    return int(mb.group(2)), md.group(1)


def collect_seeds(decker, min_overlap, max_seeds=None, span=None, span_tol=0.02):
    """Chain of slices covering as much wavelength as possible, at one slicer.

    Every usable slice of every setup is a candidate, not just the best per
    night: the 24 slices of a setup do not share a wavelength window -- the
    slicer feeds the grating at a slightly different angle per slice -- so the
    bluest and reddest slices reach further than the lowest-RMS one.
    """
    cands = []
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        parts = os.path.relpath(wc, RUNS).split(os.sep)
        night = f"{parts[0]}/{parts[2]}"
        try:
            b, dk = setup_meta(wc)
        except RuntimeError as e:
            print(f"  -- {night}: {e}")
            continue
        if dk != decker:
            print(f"  -- {night}: decker {dk}, not {decker}, skipped")
            continue
        try:
            good = rank(wc, quiet=True, span=span)
        except Exception as e:
            print(f"  -- {night}: unreadable ({type(e).__name__}), skipped")
            continue
        if not good:
            print(f"  -- {night} ({dk}): no usable slice")
            continue
        for g in good:
            g.update(file=wc, bin=b, night=night, decker=dk)
        cands.extend(good)
        print(f"  ++ {night} ({dk}): {len(good):2d} usable slices, "
              f"{min(g['wmin'] for g in good):.1f}-{max(g['wmax'] for g in good):.1f} A, "
              f"best rms {min(g['rms'] for g in good):.3f}, bin {b}")

    if not cands:
        return []

    # Widest first -- coverage is the thing that must not fall short, and RH2
    # measured that a deliberately wide seed at rms 0.296 px costs nothing a
    # throughput continuum can see.  But that argument only buys a MEANINGFUL
    # amount of coverage.  RH4's two seeds differ by 10 A in 985 (1%) while their
    # RMS differs by 0.32 vs 0.499 px, the latter sitting on pick_seed_slits'
    # 0.5 rejection boundary.  So spans within `span_tol` of the widest count as
    # equivalent and the tie goes to the lower RMS.
    widest = max(c["wmax"] - c["wmin"] for c in cands)
    tied = [c for c in cands
            if (c["wmax"] - c["wmin"]) >= widest * (1.0 - span_tol)]
    chain = [min(tied, key=lambda c: c["rms"])]
    if len(tied) > 1:
        print(f"\n{len(tied)} seeds within {100*span_tol:.0f}% of the widest "
              f"({widest:.0f} A); taking the lowest RMS of them")
    while max_seeds is None or len(chain) < max_seeds:
        end = chain[-1]["wmax"]
        reach = [c for c in cands
                 if c["wmax"] > end and c["wmin"] < end - min_overlap]
        if not reach:
            break
        chain.append(max(reach, key=lambda c: (c["wmax"], -c["rms"])))
    return chain


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--decker", default=DECKER)
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--max-seeds", type=int, default=1,
                    help="stop the chain at this many seeds.  Default 1: RH4 is "
                         "one config, and a template one exposure wide is the "
                         "shape PypeIt ships.  Width costs alignment -- "
                         "full_template cuts a window as wide as your data out "
                         "of the template at a position set by one "
                         "cross-correlation")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=(800.0, 1200.0),
                    help="physically possible wavelength span of one RH4 slice, "
                         "in A.  pick_seed_slits defaults to 400-900, which was "
                         "set from RH1 (~630 A) and RH2 (~764 A) and rejects "
                         "every RH4 slice: RH4 measures 983-994 A per exposure "
                         "at 0.476 A/px.  Span is a property of the grating, so "
                         "it has to be restated per grating")
    ap.add_argument("--span-tol", type=float, default=0.02,
                    help="seeds within this fraction of the widest span count as "
                         "equally wide, and the tie is broken on RMS")
    args = ap.parse_args()

    print(f"scanning WaveCalib files for decker={args.decker} ...")
    seeds = collect_seeds(args.decker, args.min_overlap,
                          max_seeds=args.max_seeds, span=args.span,
                          span_tol=args.span_tol)
    if not seeds:
        sys.exit("no usable RH4 seed anywhere -- nothing to build from.  Run "
                 "seed_rh4.py first, and check it used the ThArRH4 line list: "
                 "with PypeIt's KCWI default there are 4 catalogue lines in band "
                 "and nothing can solve.")

    print(f"\nstitching {len(seeds)} seed(s):")
    for s in seeds:
        print(f"   {s['night']} spat {s['spat']:4d}  {s['wmin']:.1f}-{s['wmax']:.1f} A"
              f"  rms {s['rms']:.3f}  nlines {s['nlines']}")

    cuts = []
    for a, b in zip(seeds, seeds[1:]):
        ov = a["wmax"] - b["wmin"]
        if ov <= 0:
            sys.exit(f"gap of {-ov:.0f} A between {a['night']} and {b['night']}: "
                     f"calibrate a setup with a central wavelength in between")
        cuts.append(0.5 * (b["wmin"] + a["wmax"]))
        print(f"   overlap {ov:.0f} A -> stitch at {cuts[-1]:.1f} A")

    # Written under a temporary name and moved into place: os.replace is atomic,
    # so a reduction running against the current template sees the old file or
    # the new one, never a half-written one.  The .fits extension is kept because
    # astropy infers the format from it.
    binspec = seeds[0]["bin"]
    outdir = os.path.dirname(args.out)
    tmproot = f"_tmp_{os.path.basename(args.out)}"
    tmp = os.path.join(outdir, tmproot)
    templates.build_template([s["file"] for s in seeds],
                             [s["spat"] for s in seeds],
                             cuts, binspec, tmproot, outdir=outdir,
                             ifiles=list(range(len(seeds))),
                             binning=[s["bin"] for s in seeds],
                             normalize=True, subtract_conti=True,
                             shift_wave=len(seeds) > 1, overwrite=True)

    tbl = Table.read(tmp)
    w, f = np.asarray(tbl["wave"]), np.asarray(tbl["flux"])
    d = np.diff(w)
    if not (d > 0).all():
        os.unlink(tmp)
        sys.exit("template wavelengths are not monotonic -- stitch is wrong, "
                 "the existing template is left in place")
    os.replace(tmp, args.out)
    print(f"\nwrote {args.out}")
    print(f"   {len(w)} pixels, {w.min():.1f}-{w.max():.1f} A, binspec {binspec}")
    print(f"   dispersion {np.median(d):.4f} A/px, monotonic: True")
    print(f"   flux finite: {bool(np.isfinite(f).all())}, "
          f"range {np.nanmin(f):.3g}-{np.nanmax(f):.3g}")
    print("\n   REMINDER: this template is only usable together with")
    print(f"   lamps = {os.path.join(HERE, 'ThArRH4')}")


if __name__ == "__main__":
    main()
