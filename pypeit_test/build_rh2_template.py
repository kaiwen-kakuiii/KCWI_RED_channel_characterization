#!/usr/bin/env python
"""Build the wavelength templates PypeIt does not ship for the RH2 grating.

Same hole as RH1: PypeIt 2.0.1 assigns reid_arxiv for RL, RM1, RM2 and RH3 only
(keck_kcwi.py config_specific_par), so RH1, RH2 and RH4 fall back to holy-grail,
which solves a handful of the 24 slices and leaves the flat field to die on the
rest.  See RH1_WAVELENGTH.md for the full argument; this is the RH2 instance of
it, and differs in one way that matters.

**RH2 is split across two slicers, and each needs its own template.**  The
slicer sets spectral resolution, so a template built at one slicer holds as
single features the blends another slicer resolves -- the line inventories
differ and identification fails.  That is what killed RH1's Small-slicer setup
against a Large-slicer template (0.69 A data, 1.42 A template, ratio 2.06).
RH2's nights divide cleanly:

    Large   2023-07-16 (7180), 2023-11-07 (7100), 2023-11-17 B/C (6900, 7100)
    Medium  2024-10-29, 2024-11-03, 2024-11-06, 2024-11-07 (all 6750)

so `--decker` is required rather than defaulted: building one template across
both would silently mix resolutions.

Usage:
    python build_rh2_template.py --decker Large
    python build_rh2_template.py --decker Medium --max-seeds 1
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
RUNS = os.path.join(ROOT, "RH2 pypeit run")

# One template per slicer.  Large keeps the plain name because it is the one the
# majority of the Large-slicer nights use; Medium is suffixed the way RH1's
# Small-slicer template is.
OUT = {"Large":  os.path.join(HERE, "keck_kcrm_RH2.fits"),
       "Medium": os.path.join(HERE, "keck_kcrm_RH2_medium.fits")}


def setup_meta(wavecalib_path):
    """(spectral binning, decker) of the setup a WaveCalib belongs to.

    Both come from the setup block of the .pypeit file rather than the FITS
    headers, so they describe the setup PypeIt actually grouped, not one frame.
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


def collect_seeds(decker, min_overlap, only=None, max_seeds=None, target=None):
    """Chain of slices covering as much wavelength as possible, at one slicer.

    Every usable slice of every setup is a candidate, not just the best per
    night: the 24 slices of a setup do not share a wavelength window -- the
    slicer feeds the grating at a slightly different angle per slice, spreading
    the windows over ~40 A -- so the bluest and reddest slices reach further than
    the lowest-RMS one.

    Greedy from the blue: take the bluest slice, then repeatedly take whichever
    candidate reaches furthest red while still overlapping the chain's red end by
    `min_overlap`, ties broken on RMS.  Every join is a place the wavelength
    scale can step, so the chain uses as few seeds as it can.
    """
    cands = []
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        if only and only not in wc:
            continue
        parts = os.path.relpath(wc, RUNS).split(os.sep)
        # night AND setup: two configs of one night are two arcs at different
        # central wavelengths, and labelling both by the night reads as one.
        night = f"{parts[0]}/{parts[2]}"
        try:
            b, dk = setup_meta(wc)
        except RuntimeError as e:
            print(f"  -- {night}: {e}")
            continue
        if dk != decker:
            continue
        # A reduction may be running while this rebuilds, so a WaveCalib can be
        # half-written.  Skip it rather than abort: the rest of the chain is
        # still valid and that setup can contribute on the next rebuild.
        try:
            good = rank(wc, quiet=True)
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

    if target:
        # Aim the chain at the setup being calibrated rather than starting from
        # the bluest slice: the shift full_template needs is set by how far the
        # data sits from the template's start, so aiming keeps that shift small.
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
    ap.add_argument("--decker", required=True, choices=sorted(OUT),
                    help="slicer to build for; required because the slicer sets "
                         "spectral resolution and a template must match it")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--only", default=None,
                    help="substring a WaveCalib path must contain, e.g. a single "
                         "setup -- with --max-seeds 1 this builds a template one "
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
    out = args.out or OUT[args.decker]

    print(f"scanning WaveCalib files for decker={args.decker} ...")
    seeds = collect_seeds(args.decker, args.min_overlap, only=args.only,
                          max_seeds=args.max_seeds, target=args.target)
    if not seeds:
        sys.exit(f"no usable {args.decker} seed anywhere -- nothing to build from")

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
    # may be running against the current template while this rebuilds it:
    # os.replace is atomic, so a concurrent reader sees the old file or the new
    # one, never the half-written one.  The .fits extension is kept because
    # astropy infers the format from it.
    binspec = seeds[0]["bin"]
    outdir = os.path.dirname(out)
    tmproot = f"_tmp_{os.path.basename(out)}"
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
    os.replace(tmp, out)
    print(f"\nwrote {out}")
    print(f"   {len(w)} pixels, {w.min():.1f}-{w.max():.1f} A, binspec {binspec}")
    print(f"   dispersion {np.median(d):.4f} A/px, monotonic: True")
    print(f"   flux finite: {bool(np.isfinite(f).all())}, "
          f"range {np.nanmin(f):.3g}-{np.nanmax(f):.3g}")


if __name__ == "__main__":
    main()
