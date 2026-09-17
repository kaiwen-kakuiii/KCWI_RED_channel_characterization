#!/usr/bin/env python
"""Build wavelength templates for the RM1 central wavelengths PypeIt's own does not reach.

PypeIt ships `keck_kcrm_RM1.fits` and lists RM1 as a supported setup.  Measured
here, that template is real and usable -- but only over part of this run:

    shipped keck_kcrm_RM1.fits    6293.1 - 7733.7 A, 0.7006 A/px, centre 7013
    an RM1 exposure               1446 A wide (2064 binned px x 0.7006)

so the fraction of detector with template under it runs from 74% at cenwave
7390 down to 39% at 6130.  Gated on three configs, the failure tracks that
number and reproduces RH3's ~50% threshold exactly:

    2023-12-10 C   cenwave 6630   73% overlap   24/24 solved   rms median 0.200 px
    2024-03-15 B   cenwave 6300   50% overlap   13/24          rms median 0.982 px
    2024-05-02 B   cenwave 6200   44% overlap   15/24          rms median 0.218 px

**Neither broken config aborted.**  PypeIt logged no `Not enough useful IDs` and
wrote a WaveCalib for all 24 slices either way; the bad ones simply carry a
dispersion of 0.48-0.91 A/px against the 0.713 the good slices agree on.  A
template covering half your detector does not half-work, and it does not say so.

## What splits a template here

**Central wavelength**, as on RH3 -- seven of them, 6130 to 7510, which is 1380 A
of grating rotation.  **And slicer**, as on RH2: the slicer sets spectral
resolution, and a template built at one resolution holds as single features the
blends another resolves.  RM1 is the first run in this project where BOTH split
at once -- five nights Large/2x2, one night (2024-04-01) Small/1x1 -- so
`collect_seeds` refuses to mix deckers rather than warning about it.

## One template serves four central wavelengths

RH3 built one template per cenwave because it had two.  Seven would be absurd,
and unnecessary: an RM1 exposure is 1446 A wide, so a template seeded at 6300
covers 5577-7023 A and sits under

    6130   88%       6200   92%       6300  100%       6480   88%

of each detector -- far above the ~50% where the shipped one fails.  So the blue
group takes ONE template.  6630, 7390 and 7510 keep the shipped template, which
gates clean there.

## Where the seeds come from

The RH3 section-4 bootstrap.  A config the shipped template FAILS to calibrate
still produces slices that pass the strict seed criteria -- 6 of 24 on
2024-03-15, 10 of 24 on 2024-05-02, dispersions agreeing to a fraction of a
percent and blue ends where a 6300 setup must sit.  Those seed the new template,
and the config then gates against it.  "FAIL 13/24" reads like a dead end and is
not one.

## --span is not optional here

`pick_seed_slits.rank()` defaults to a 400-900 A window for one slice, set from
RH1 (~630 A) and RH2 (~764 A).  An RM1 slice spans **~1453 A**.  Called with the
default every single slice is rejected as nonsense and the build reports "no
usable seed", which looks like bad data and is a unit mismatch.  This script
therefore defaults `--span` to 1200-1700 and passes it through.

Usage:
    python build_rm1_template.py --cenwave 6300
    python build_rm1_template.py --cenwave 6130 --decker Small
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
RUNS = os.path.join(ROOT, "RM1 pypeit run")

# Measured: RM1 runs ~0.7060 A/px over 2064 binned px = ~1457 A per SLICE at
# 2x2, centred near RCWAVE.  The 24 slices together span a little more, because
# each is fed the grating at a slightly different angle -- on 2023-12-10 C the
# blue ends run 5852.5-5991.9 A, a 139 A union that alternates odd/even slice,
# which is also why RM1 needs a shift_tol above the 50 A default (see
# gate_rm1_template.py).
HALF_SPAN = 728.0
SPAN_DEFAULT = (1200.0, 1700.0)


def outname(cenwave, decker):
    """Never `keck_kcrm_RM1.fits` -- that name belongs to PypeIt's shipped file.

    The runs point reid_arxiv at an absolute path in this repo, so a clash would
    not silently resolve to the wrong file; a reader should not have to know
    that to tell ours from PypeIt's.  Decker is in the name only when it is not
    the run's majority Large, so the common case reads simply.
    """
    tail = "" if decker.lower() == "large" else f"_{decker.lower()}"
    return os.path.join(HERE, f"keck_kcrm_RM1_{cenwave}{tail}.fits")


def setup_meta(wavecalib_path):
    """(spectral binning, decker, cenwave) of the setup a WaveCalib belongs to.

    All three come from the setup block of the .pypeit file rather than from a
    frame header, so they describe the setup PypeIt actually grouped.  Cenwave
    is rounded to 10 A, not 50 as RH3 used: RM1's settings sit 100-130 A apart
    (6130, 6200, 6300), and a 50 A round merges 6130 with nothing but would
    misreport the spacing this run turns on.
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
    # KCWI writes binning as spatial,spectral.
    return int(mb.group(2)), md.group(1), int(round(float(mc.group(1)) / 10.0) * 10)


def collect_seeds(cenwaves, decker, min_overlap, span, only=None, max_seeds=None):
    """Chain of slices covering as much wavelength as possible.

    Every usable slice of every matching setup is a candidate, not just the best
    per night: the 24 slices of a setup do not share a wavelength window, so the
    bluest and reddest slices reach further than the lowest-RMS one.

    Greedy from the bluest seed, then repeatedly take whichever candidate reaches
    furthest red while still overlapping the chain's red end by `min_overlap`,
    ties broken on RMS.  Every join is a place the wavelength scale can step, so
    the chain uses as few seeds as it can.

    `cenwaves` is a LIST, unlike RH3's single value.  RM1's blue group sits at
    6130/6200/6300/6480 within one 1446 A exposure width, so slices from
    neighbouring settings are legitimate seeds for one another and using them
    widens the chain where a single setting would leave it short.
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
        if cw not in cenwaves or dk != decker:
            continue
        deckers.add(dk)
        try:
            good = rank(wc, quiet=True, span=span)
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
    # other slicer resolves.  RM1 has both, so this guard can actually fire.
    if len(deckers) > 1:
        sys.exit(f"seeds span deckers {sorted(deckers)} -- a single template "
                 f"cannot serve two spectral resolutions (RH2.md)")
    if not cands:
        return []

    # Start from the BLUEST seed.  RH3 learned this the hard way: the 24 slices
    # each record a different window, so a template aimed mid-range leaves the
    # blue slices with nothing under them and they go unsolved.
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
    ap.add_argument("--cenwave", required=True, type=int,
                    help="central wavelength this template is FOR; also its name")
    ap.add_argument("--seed-cenwaves", type=int, nargs="+", default=None,
                    help="which settings may contribute seeds.  Defaults to the "
                         "target alone.  RM1's blue group (6130/6200/6300/6480) "
                         "sits inside one 1446 A exposure width, so neighbours "
                         "are legitimate seeds and widen the chain")
    ap.add_argument("--decker", default="Large",
                    help="slicer.  A template cannot serve two: the slicer sets "
                         "spectral resolution, and one built at Large holds as "
                         "single features the blends Small resolves (RH2.md)")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=list(SPAN_DEFAULT),
                    help="physically possible span of ONE RM1 slice, in A.  "
                         "Measured ~1453 A.  The pick_seed_slits default is "
                         "400-900 A, set from RH1/RH2, and with it every RM1 "
                         "slice is rejected as nonsense")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--only", default=None,
                    help="substring a WaveCalib path must contain, e.g. one setup")
    ap.add_argument("--max-seeds", type=int, default=None,
                    help="stop the chain at this many seeds.  Width costs "
                         "alignment: full_template cuts a window as wide as your "
                         "data at a position set by ONE cross-correlation, so a "
                         "template much wider than an exposure offers several "
                         "plausible alignments of the same line forest.  That "
                         "cost RH1 four slices")
    args = ap.parse_args()

    out = args.out or outname(args.cenwave, args.decker)
    span = (float(args.span[0]), float(args.span[1]))
    seed_cw = args.seed_cenwaves or [args.cenwave]

    print(f"scanning RM1 WaveCalib files: decker {args.decker}, seed cenwaves "
          f"{seed_cw}, slice span window {span[0]:.0f}-{span[1]:.0f} A ...")
    seeds = collect_seeds(seed_cw, args.decker, args.min_overlap, span,
                          only=args.only, max_seeds=args.max_seeds)
    if not seeds:
        sys.exit(f"no usable seed for cenwave {args.cenwave} ({args.decker}) -- "
                 f"nothing to build from")

    print(f"\nstitching {len(seeds)} seed(s):")
    for s in seeds:
        print(f"   {s['night']} (cw {s['cenwave']}) spat {s['spat']:4d}  "
              f"{s['wmin']:.1f}-{s['wmax']:.1f} A  rms {s['rms']:.3f}  "
              f"disp {s['disp']:.4f}")

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

    # Provenance in the FILE, not inferred later.  gate_rm1_template.py reads
    # these to say whether a gate is generalisation or reproduction; RH3 section
    # 9 shows what happens when it has to infer -- once a night is gated its
    # WaveCalib passes rank() too, so every night looks like a seed source and
    # a reproduction gets reported as independent.
    from astropy.io import fits as _fits
    with _fits.open(tmp) as hdu:
        hdr = hdu[0].header
        hdr["SEEDNITE"] = (",".join(sorted({s["night"] for s in seeds}))[:68],
                           "setups whose slices seeded this template")
        hdr["SEEDSPAT"] = (",".join(str(s["spat"]) for s in seeds)[:68],
                           "spat_ids of the seed slices")
        hdr["NSEED"] = (len(seeds), "number of seed slices stitched")
        hdr["CENWAVE"] = (args.cenwave, "central wavelength this template is for")
        hdr["DECKER"] = (args.decker, "slicer; a template cannot serve two")
        hdu.writeto(tmp, overwrite=True)

    os.replace(tmp, out)
    d = Table.read(out)
    w = np.asarray(d["wave"], float)
    print(f"\nwrote {out}\n   {w.min():.1f}-{w.max():.1f} A, {w.size} px, "
          f"{np.median(np.diff(w)):.4f} A/px, binspec {binspec}")
    print(f"   an RM1 exposure is {2*HALF_SPAN:.0f} A wide, so this covers "
          f"{100*min(1.0,(w.max()-w.min())/(2*HALF_SPAN)):.0f}% of one")


if __name__ == "__main__":
    main()
