#!/usr/bin/env python
"""Build wavelength templates for the RM2 central wavelengths PypeIt's own does not reach.

PypeIt ships `keck_kcrm_RM2.fits` and lists RM2 as a supported setup.  Measured
here, that template is real and usable -- at ONE of the five central wavelengths
in this run:

    shipped keck_kcrm_RM2.fits    8608.0 - 10595.0 A, 0.9659 A/px, centre 9601
    an RM2 exposure               1994 A wide (2064 binned px x 0.9659)

    cenwave 9850   Medium 2x2   87% overlap   24/24 solved, disp agree 0.4%   CLEAN
    cenwave 8950   Medium 2x2   67%            2/24 solved, flat CRASHED
    cenwave 8900   Medium 2x2   65%
    cenwave 8850   Small  1x1   62%
    cenwave 7750   Small  1x1    7%

RM1 measured 65% not enough and 74% not enough across a slicer change; RM2 puts
the Medium threshold between 67% and 87% and shows the 67% failure is not quiet
the way RM1's was -- at 8950 `full_template` logged `Not enough useful IDs` on 22
of 24 slices and the flat field died with `ValueError: zero-size array`, which is
RH1's accidental safety net firing.

## The seeds are richer than the setup that passed

The two slices that DID solve at 8950 fit **25-28 lines** at rms 0.030-0.057 over
8016-10067 A, against **10-14 lines** on every slice of the clean 9850 setup.
That is not a paradox: the FeAr arc against PypeIt's default FeI/ArI/ArII
catalogue thins out toward the red, so a bluer detector sees more of it.  It is
also why `--nline-min` defaults to 8 here rather than pick_seed_slits' 20 -- at
20, the setup this project would most like to seed from reports zero usable
slices, and the bootstrap cannot start.  What carries the weight at these line
counts is dispersion agreement across the 24 slices and between independent
nights, never the rms (RH1_PROCEDURE.md's "few points, small residual").

## One Medium template serves 8900 and 8950

Seeded from 2024-06-11 B's two good slices, a template covering 8016-10022 A sits
under ~99% of an 8950 detector (7953-9947) and ~94% of an 8900 one (7903-9897).
Both 8900 nights and the other two 8950 nights are then INDEPENDENT tests of it,
which is the property RH3 section 9 says to preserve deliberately: seed from one
night, gate on the others.

9850 keeps PypeIt's own template, which gates clean there.

## The Small configs need their own, as always

7750 and 8850 are Small/1x1.  A template cannot serve two slicers -- the slicer
sets spectral resolution, and one built at Medium holds as single features the
blends Small resolves (RH2.md) -- so `collect_seeds` refuses to mix deckers
rather than warning.  7750 is the hard case at 7% shipped coverage: 6753-8747 A
is entirely outside the shipped template and mostly outside any Medium seed we
have, so it must bootstrap from its own slices or from holy-grail.

## --span is not optional here

`pick_seed_slits.rank()` defaults to a 400-900 A window for one slice (RH1 ~630,
RH2 ~764); RM1 had to widen it to 1200-1700 for its ~1453 A slices.  An RM2 slice
spans **~1994 A** and is rejected as nonsense by both.  Default here is
1700-2300, and it is passed through.

Usage:
    python build_rm2_template.py --cenwave 8900 --seed-cenwaves 8950
    python build_rm2_template.py --cenwave 7750 --decker Small
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
RUNS = os.path.join(ROOT, "RM2 pypeit run")

# Measured: RM2 runs 0.9585-0.9764 A/px over 2064 binned px = ~1970-2010 A per
# SLICE at 2x2, centred near RCWAVE.  The 24 slices together span a little more,
# because each is fed the grating at a slightly different angle -- on 2023-09-23 B
# the blue ends run 8822.5-9029.2 A, a 207 A union that alternates odd/even slice,
# which is why RM2 needs a shift_tol of 150 A (see gate_rm2_template.py).
HALF_SPAN = 997.0
SPAN_DEFAULT = (1700.0, 2300.0)

# Lines fitted per slice, and this is a property of the GRATING and the LAMP.
# pick_seed_slits' 20 comes from RM1/RH, where FeAr fits 70-80 lines in a slice;
# RM2's 9850 setup fits 10-14 and its 8950 seeds 25-28, so 20 would reject the
# clean setup outright.  See the note in pick_seed_slits.py.
NLINE_MIN_DEFAULT = 8


def outname(cenwave, decker):
    """Never `keck_kcrm_RM2.fits` -- that name belongs to PypeIt's shipped file.

    The runs point reid_arxiv at an absolute path in this repo, so a clash would
    not silently resolve to the wrong file; a reader should not have to know
    that to tell ours from PypeIt's.  Decker is in the name only when it is not
    the run's majority Medium, so the common case reads simply.
    """
    tail = "" if decker.lower() == "medium" else f"_{decker.lower()}"
    return os.path.join(HERE, f"keck_kcrm_RM2_{cenwave}{tail}.fits")


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


def collect_seeds(cenwaves, decker, min_overlap, span, only=None, max_seeds=None,
                  nline_min=NLINE_MIN_DEFAULT):
    """Chain of slices covering as much wavelength as possible.

    Every usable slice of every matching setup is a candidate, not just the best
    per night: the 24 slices of a setup do not share a wavelength window, so the
    bluest and reddest slices reach further than the lowest-RMS one.

    Greedy from the bluest seed, then repeatedly take whichever candidate reaches
    furthest red while still overlapping the chain's red end by `min_overlap`,
    ties broken on RMS.  Every join is a place the wavelength scale can step, so
    the chain uses as few seeds as it can.

    `cenwaves` is a LIST, unlike RH3's single value.  RM2's 8850/8900/8950 sit
    100 A apart inside one 1994 A exposure width, so slices from neighbouring
    settings are legitimate seeds for one another -- which is the whole point
    here, since 8950 is the only setting whose own slices solved.
    """
    cands, deckers = [], set()
    for wc in sorted(glob.glob(os.path.join(RUNS, "*", "pypeit_run", "keck_kcrm_*",
                                            "Calibrations", "WaveCalib_*.fits"))):
        # `only` is a LIST here, unlike build_rm1_template.py's single string.
        # RM2 needs it that way: the Medium template is seeded from 2024-01-04
        # (the only slice reaching 7871 A) plus 2024-06-11, and 2025-01-01 /
        # 2025-01-02 / 2024-04-30 are deliberately WITHHELD so that gating them
        # is an independent test rather than a reproduction (RH3 section 9).
        # With a single pattern the withheld nights could not be expressed.
        if only and not any(o in wc for o in only):
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
            good = rank(wc, quiet=True, span=span, nline_min=nline_min)
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
                         "target alone.  RM2's 8850/8900/8950 sit inside one "
                         "1994 A exposure width, so neighbours are legitimate "
                         "seeds -- and at 8900 they are the ONLY ones, since "
                         "only 8950 produced solved slices")
    ap.add_argument("--decker", default="Medium",
                    help="slicer.  A template cannot serve two: the slicer sets "
                         "spectral resolution, and one built at Medium holds as "
                         "single features the blends Small resolves (RH2.md)")
    ap.add_argument("--span", type=float, nargs=2, metavar=("MIN", "MAX"),
                    default=list(SPAN_DEFAULT),
                    help="physically possible span of ONE RM2 slice, in A.  "
                         "Measured 1958-2011 A.  The pick_seed_slits default is "
                         "400-900 A and RM1's is 1200-1700; with either, every "
                         "RM2 slice is rejected as nonsense")
    ap.add_argument("--nline-min", type=int, default=NLINE_MIN_DEFAULT,
                    help="lines a slice must fit to seed a template.  8 here "
                         "against pick_seed_slits' 20: RM2's FeAr arc fits "
                         "10-14 lines per slice at 9850 and 25-28 at 8950, so "
                         "20 rejects the setup that gated 24/24")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-overlap", type=float, default=100.0,
                    help="A a seed must overlap the chain before it")
    ap.add_argument("--only", nargs="+", default=None,
                    help="substrings a WaveCalib path may contain; a path "
                         "matching ANY of them is a seed candidate.  Use it to "
                         "WITHHOLD nights so gating them stays an independent "
                         "test of the template rather than a reproduction")
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

    print(f"scanning RM2 WaveCalib files: decker {args.decker}, seed cenwaves "
          f"{seed_cw}, slice span window {span[0]:.0f}-{span[1]:.0f} A, "
          f"nline_min {args.nline_min} ...")
    seeds = collect_seeds(seed_cw, args.decker, args.min_overlap, span,
                          only=args.only, max_seeds=args.max_seeds,
                          nline_min=args.nline_min)
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

    # Provenance in the FILE, not inferred later.  gate_rm2_template.py reads
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
    print(f"   an RM2 exposure is {2*HALF_SPAN:.0f} A wide, so this covers "
          f"{100*min(1.0,(w.max()-w.min())/(2*HALF_SPAN)):.0f}% of one")


if __name__ == "__main__":
    main()
