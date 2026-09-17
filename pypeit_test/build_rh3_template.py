#!/usr/bin/env python
"""Build wavelength templates for RH3, whose shipped one is at the wrong cenwave.

RH1, RH2 and RH4 had no `reid_arxiv` at all.  **RH3 does** -- `keck_kcwi.py
config_specific_par` maps dispname RH3 to `keck_kcrm_RH3.fits`, and PypeIt lists
RH3 as a supported setup.  That template is real and presumably correct for the
data it was built from; it simply does not cover this data:

    shipped keck_kcrm_RH3.fits    8599.4 - 9389.0 A
    these nights, cenwave 8600      ~8204 - 8996 A     ~50% overlap
    these nights, cenwave 8400      ~8004 - 8796 A     ~25% overlap

Measured consequence (seed_rh3.py's docstring has the log): 1 usable slice of 24
on the 8600 night and 10 on the 8400 night, but every one of them at rms 0.000
through too few points, so **0 slices meet the seed criteria on either**.  A
template that covers half your detector does not half-work; it fails.

**The split here is CENTRAL WAVELENGTH, not slicer.**  RH2 needed one template
per slicer because the slicer sets spectral resolution and a template built at
one resolution holds as single features the blends another resolves.  Every RH3
night is Medium/2x2, so resolution is constant and that argument does not apply.
What varies is cenwave -- 8600 on three nights, 8400 on one -- and two setups
200 A apart do not share a line inventory at the ends.  Hence one template per
cenwave, aimed with `--target`, each about one exposure wide.

**One exposure wide is deliberate**, and it is RH4's "good case": `full_template`
cuts a window as wide as your data out of the template at a position set by a
single cross-correlation, so a template much wider than one exposure offers
several plausible alignments of the same line forest and can be cut in the wrong
place.  That ambiguity cost RH1 four slices (RH1_PROCEDURE.md problem 3).  A
chain across both cenwaves would be ~990 A against a 793 A exposure and would
reintroduce exactly that, to no benefit -- each setup only ever needs its own
793 A.

The output names carry the cenwave and **never collide with the shipped
`keck_kcrm_RH3.fits`**.  The runs point `reid_arxiv` at an absolute path in this
repo, so a name clash would not silently resolve to the wrong file, but a reader
should not have to know that to tell ours from PypeIt's.

Usage:
    python build_rh3_template.py --cenwave 8600
    python build_rh3_template.py --cenwave 8400 --max-seeds 1
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
RUNS = os.path.join(ROOT, "RH3 pypeit run")

# One template per central wavelength.  Both are suffixed: the unsuffixed name
# belongs to PypeIt's shipped template and reusing it would make two different
# files indistinguishable in a log line.
OUT = {8600: os.path.join(HERE, "keck_kcrm_RH3_8600.fits"),
       8400: os.path.join(HERE, "keck_kcrm_RH3_8400.fits")}

# Measured detector coverage: RH3 runs 0.4045 A/px over 2064 binned px = 835 A
# per SLICE, centred on RCWAVE.  The 24 slices together span ~897 A, because each
# is fed the grating at a slightly different angle.  (The shipped
# keck_kcrm_RH3.fits implies 0.384 A/px, which is one more sign it was built from
# a different setup than these nights.)
HALF_SPAN = 417.0


def setup_meta(wavecalib_path):
    """(spectral binning, decker, cenwave) of the setup a WaveCalib belongs to.

    All three come from the setup block of the .pypeit file rather than from a
    frame header, so they describe the setup PypeIt actually grouped.
    """
    sdir = os.path.dirname(os.path.dirname(wavecalib_path))
    pfile = glob.glob(os.path.join(sdir, "*.pypeit"))
    if not pfile:
        raise RuntimeError(f"no .pypeit beside {wavecalib_path}")
    txt = open(pfile[0]).read()
    mb = re.search(r"binning:\s*(\d+),(\d+)", txt)
    md = re.search(r"decker:\s*(\S+)", txt)
    mc = re.search(r"cenwave:\s*([0-9.]+)", txt)
    if not (mb and md and mc):
        raise RuntimeError(f"cannot read binning/decker/cenwave for {wavecalib_path}")
    # KCWI writes binning as spatial,spectral.  Round cenwave: the same nominal
    # setting is written 8599.9794922 and 8599.98 by different frames.
    return int(mb.group(2)), md.group(1), int(round(float(mc.group(1)) / 50.0) * 50)


def collect_seeds(cenwave, min_overlap, only=None, max_seeds=None, target=None):
    """Chain of slices covering as much wavelength as possible, at one cenwave.

    Every usable slice of every matching setup is a candidate, not just the best
    per night: the 24 slices of a setup do not share a wavelength window -- the
    slicer feeds the grating at a slightly different angle per slice -- so the
    bluest and reddest slices reach further than the lowest-RMS one.

    Greedy from the seed that overlaps `target` most, then repeatedly take
    whichever candidate reaches furthest red while still overlapping the chain's
    red end by `min_overlap`, ties broken on RMS.  Every join is a place the
    wavelength scale can step, so the chain uses as few seeds as it can.
    """
    cands, deckers = [], set()
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        if only and only not in wc:
            continue
        parts = os.path.relpath(wc, RUNS).split(os.sep)
        night = f"{parts[0]}/{parts[2]}"
        try:
            b, dk, cw = setup_meta(wc)
        except RuntimeError as e:
            print(f"  -- {night}: {e}")
            continue
        if cw != cenwave:
            continue
        deckers.add(dk)
        try:
            good = rank(wc, quiet=True)
        except Exception as e:
            print(f"  -- {night}: unreadable ({type(e).__name__}), skipped")
            continue
        if not good:
            print(f"  -- {night} ({dk}, {cw}): no usable slice")
            continue
        for g in good:
            g.update(file=wc, bin=b, night=night, decker=dk, cenwave=cw)
        cands.extend(good)
        print(f"  ++ {night} ({dk}, {cw}): {len(good):2d} usable slices, "
              f"{min(g['wmin'] for g in good):.1f}-{max(g['wmax'] for g in good):.1f} A, "
              f"best rms {min(g['rms'] for g in good):.3f}, bin {b}")

    # The RH2 lesson, enforced rather than assumed: mixing slicers mixes spectral
    # resolutions, and the template then holds as single features the blends the
    # other slicer resolves.  RH3 should be all Medium; say so if it ever is not.
    if len(deckers) > 1:
        sys.exit(f"cenwave {cenwave} spans deckers {sorted(deckers)} -- a single "
                 f"template cannot serve two spectral resolutions (RH2.md); "
                 f"split on decker as build_rh2_template.py does")
    if not cands:
        return []

    # Start from the BLUEST seed unless a target says otherwise.
    #
    # This is where RH3 parts company with RH4's "one exposure wide is the good
    # case".  That holds when the template only has to cover one slice.  It does
    # not here: the 24 slices of a setup each record a different window -- RH3's
    # run 8136.4-9033.0 A, a 897 A union -- while any ONE slice spans only 832 A.
    # A single-seed template is therefore 65 A short of its own setup no matter
    # which seed is chosen, and the slices at whichever end it misses go
    # unsolved.  The chain has to walk blue to red across the union, which is
    # what a blue start does; aiming it with --target starts it mid-range and
    # silently leaves the blue slices uncovered.
    if target is None:
        chain = [min(cands, key=lambda c: (c["wmin"], c["rms"]))]
    else:
        lo, hi = target
        chain = [max(cands, key=lambda c: (min(c["wmax"], hi) - max(c["wmin"], lo),
                                           -c["rms"]))]
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
    ap.add_argument("--cenwave", required=True, type=int, choices=sorted(OUT),
                    help="central wavelength to build for; required because two "
                         "setups 200 A apart do not share a line inventory at "
                         "the ends")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--only", default=None,
                    help="substring a WaveCalib path must contain, e.g. a single "
                         "setup -- with --max-seeds 1 this builds a template one "
                         "exposure wide, the shape PypeIt ships")
    ap.add_argument("--target", type=float, nargs=2, metavar=("LO", "HI"),
                    default=None,
                    help="wavelength range to aim the chain at.  OMIT IT for "
                         "normal use: the chain then starts at the bluest seed "
                         "and walks red across the whole union of slice windows, "
                         "which is what a setup needs.  Aiming starts the chain "
                         "mid-range and leaves the bluest slices uncovered")
    ap.add_argument("--max-seeds", type=int, default=None,
                    help="stop the chain at this many seeds.  Width costs "
                         "alignment -- see the module docstring")
    args = ap.parse_args()
    out = args.out or OUT[args.cenwave]
    target = tuple(args.target) if args.target else None

    print(f"scanning WaveCalib files for cenwave={args.cenwave} "
          + (f"(aimed at {target[0]:.0f}-{target[1]:.0f} A) ..." if target
             else "(blue-start chain across all slice windows) ..."))
    seeds = collect_seeds(args.cenwave, args.min_overlap, only=args.only,
                          max_seeds=args.max_seeds, target=target)
    if not seeds:
        sys.exit(f"no usable seed at cenwave {args.cenwave} -- nothing to build from")

    print(f"\nstitching {len(seeds)} seed(s):")
    for s in seeds:
        print(f"   {s['night']} spat {s['spat']:4d}  {s['wmin']:.1f}-{s['wmax']:.1f} A"
              f"  rms {s['rms']:.3f}")

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
    # one, never the half-written one.
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

    # Stamp the provenance INTO the file.  Without it, "was this night one of the
    # seeds?" has to be answered by looking at which WaveCalibs currently pass
    # rank() -- and after a night is gated, its WaveCalib passes too, so every
    # night looks like a seed source and a gate cannot tell reproduction from
    # generalisation.  That misreported RH3's own 2024-10-29 gate as independent
    # when the template was built from exactly that night.
    from astropy.io import fits as _fits
    with _fits.open(tmp, mode="update") as _h:
        hd = _h[1].header if len(_h) > 1 else _h[0].header
        hd["SEEDNITE"] = (",".join(sorted({s["night"].split("/")[0]
                                           for s in seeds})),
                          "nights whose slices seeded this template")
        hd["SEEDSPAT"] = (",".join(str(s["spat"]) for s in seeds),
                          "spat_ids of the seed slices")
        hd["NSEED"] = (len(seeds), "number of stitched seeds")
        hd["CENWAVE"] = (args.cenwave, "central wavelength this template serves")
    os.replace(tmp, out)
    print(f"\nwrote {out}")
    print(f"   {len(w)} pixels, {w.min():.1f}-{w.max():.1f} A, binspec {binspec}")
    print(f"   dispersion {np.median(d):.4f} A/px, monotonic: True")
    print(f"   flux finite: {bool(np.isfinite(f).all())}, "
          f"range {np.nanmin(f):.3g}-{np.nanmax(f):.3g}")


if __name__ == "__main__":
    main()
