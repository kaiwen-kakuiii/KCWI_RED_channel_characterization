#!/usr/bin/env python
"""Build the wavelength template PypeIt does not ship for the RH1 grating.

PypeIt 2.0.1 has reid_arxiv templates for RL, RM1, RM2 and RH3 only
(keck_kcwi.py config_specific_par).  For RH1 it falls back to holy-grail, which
here solved 1 of 24 slices per setup and left the rest unusable -- and an
unsolved slice kills the IFU flat field (flatfield.py illum_profile_spectral
takes np.min of an empty wavelength array).

One RH1 exposure covers only ~630 A, so no single setup spans the 5824-7000 A
these nights need.  The template is therefore stitched from the best slice of
several setups at different central wavelengths, which is what
pypeit.core.wavecal.templates.build_template exists for.

Seeds are chosen by pick_seed_slits.py, one per distinct central wavelength, and
stitched at the midpoint of each consecutive pair's overlap.

Usage:
    python build_rh1_template.py [--out <path>] [--min-overlap 50]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_seed_slits import rank                       # noqa: E402
from pypeit.core.wavecal import templates              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "RH1 pypeit run")


def setup_binning(wavecalib_path):
    """Spectral binning of the setup a WaveCalib belongs to, from its .pypeit file."""
    sdir = os.path.dirname(os.path.dirname(wavecalib_path))
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if pfile:
        m = re.search(r"binning:\s*(\d+),(\d+)", open(pfile[0]).read())
        if m:
            # KCWI writes binning as spatial,spectral
            return int(m.group(2))
    raise RuntimeError(f"cannot determine binning for {wavecalib_path}")


def collect_seeds(min_overlap, only=None, max_seeds=None, target=None):
    """Choose a chain of slices covering as much wavelength as possible.

    Every usable slice of every setup is a candidate, not just the best one per
    night: the 24 slices of one setup do not share a wavelength window -- the
    slicer feeds the grating at a slightly different angle per slice, spreading
    their windows over ~40 A -- so the bluest and reddest slices reach further
    than the one with the lowest RMS.

    Selection is greedy from the blue: take the bluest slice, then repeatedly
    take whichever candidate reaches furthest to the red while still overlapping
    the chain's current red end by `min_overlap`, breaking ties on RMS.  Each
    stitch therefore has real overlap to align on, and the chain uses as few
    seeds as possible -- every join is a place the wavelength scale can step.
    """
    cands = []
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        if only and only not in wc:
            continue
        # night AND setup: two configs of the same night are two different arcs
        # at different central wavelengths, and labelling both by the night alone
        # makes a two-frame stitch read as one frame.
        parts = os.path.relpath(wc, RUNS).split(os.sep)
        night = f"{parts[0]}/{parts[2]}"
        # A reduction may be running while this rebuilds, so a WaveCalib can be
        # half-written.  Skip it rather than abort: the rest of the chain is
        # still valid, and that setup can contribute on the next rebuild.
        try:
            good = rank(wc, quiet=True)
        except Exception as e:
            print(f"  -- {night}: unreadable ({type(e).__name__}), skipped")
            continue
        if not good:
            print(f"  -- {night}: no usable slice")
            continue
        b = setup_binning(wc)
        for g in good:
            g.update(file=wc, bin=b, night=night)
        cands.extend(good)
        print(f"  ++ {night}: {len(good):2d} usable slices, "
              f"{min(g['wmin'] for g in good):.1f}-{max(g['wmax'] for g in good):.1f} A, "
              f"best rms {min(g['rms'] for g in good):.3f}, bin {b}")

    if not cands:
        return []

    if target:
        # Aim the chain at the setup being calibrated rather than starting from
        # the bluest slice available: the shift full_template needs is set by how
        # far the data sits from the template's start, so a template aimed at the
        # data keeps that shift small.
        lo, hi = target
        chain = [max(cands, key=lambda c: (min(c["wmax"], hi) - max(c["wmin"], lo),
                                           -c["rms"]))]
    else:
        chain = [min(cands, key=lambda c: (c["wmin"], c["rms"]))]
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
    ap.add_argument("--out", default=os.path.join(ROOT, "pypeit_test", "keck_kcrm_RH1.fits"))
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--only", default=None,
                    help="substring a WaveCalib path must contain, e.g. a single "
                         "setup -- use with --max-seeds 1 to build a template one "
                         "exposure wide, the shape PypeIt ships")
    ap.add_argument("--target", type=float, nargs=2, metavar=("LO", "HI"),
                    default=None,
                    help="wavelength range of the setup to be calibrated; the "
                         "chain starts from the seed overlapping it most")
    ap.add_argument("--max-seeds", type=int, default=None,
                    help="stop the chain at this many seeds.  Width costs "
                         "alignment: full_template cuts a window as wide as your "
                         "data out of the template at a position set by one "
                         "cross-correlation, so a template much wider than one "
                         "exposure offers many plausible alignments of the same "
                         "line forest and can be cut in the wrong place")
    args = ap.parse_args()

    print("scanning WaveCalib files ...")
    seeds = collect_seeds(args.min_overlap, only=args.only,
                          max_seeds=args.max_seeds, target=args.target)
    if not seeds:
        sys.exit("no usable seed anywhere -- nothing to build from")

    print(f"\nstitching {len(seeds)} seed(s):")
    for s in seeds:
        print(f"   {s['night']} spat {s['spat']:4d}  {s['wmin']:.1f}-{s['wmax']:.1f} A")

    cuts = []
    for a, b in zip(seeds, seeds[1:]):
        ov = a["wmax"] - b["wmin"]
        if ov <= 0:
            sys.exit(f"gap of {-ov:.0f} A between {a['night']} and {b['night']}: "
                     f"calibrate a setup with a central wavelength in between")
        cuts.append(0.5 * (b["wmin"] + a["wmax"]))
        print(f"   overlap {ov:.0f} A -> stitch at {cuts[-1]:.1f} A")

    # Written under a temporary name and moved into place, because a reduction
    # may be running against the current template while this rebuilds it: os.replace
    # is atomic, so a concurrent reader sees either the old file or the new one,
    # never the half-written one.
    binspec = seeds[0]["bin"]
    outdir = os.path.dirname(args.out)
    # Keeps the .fits extension: astropy infers the format from it, and a name
    # ending .fits.new is not a FITS file as far as Table.read is concerned.
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


if __name__ == "__main__":
    main()
