#!/usr/bin/env python
"""Rank the slits of a WaveCalib by how trustworthy their solution is.

Holy-grail (no reid_arxiv template) solves only a handful of the 24 KCRM slices
and fails the rest silently -- a "fit" with 6 lines and a 10^5 A span is a
returned object, not a solution.  This picks the ones fit to seed a template.

A slit passes only if all of:
    rms      < 0.5 binned px      (the fit itself)
    nlines  >= 20                 (enough constraints across the range)
    span     = 400-900 A          (RH1 covers ~630 A; anything else is nonsense)
    monotonic dispersion, and d(lambda)/dpix within 25% of the median

Usage:  python pick_seed_slits.py [--span MIN MAX] [--shift-tol A]
                                 [--nline-min N] <WaveCalib_*.fits> [...]
"""
import glob
import os
import sys

import numpy as np
from pypeit.slittrace import SlitTraceSet
from pypeit.wavecalib import WaveCalib

RMS_MAX, NLINE_MIN, SPAN_MIN, SPAN_MAX, DISP_TOL = 0.5, 20, 400.0, 900.0, 0.25

# How far a slice's blue end may sit from the median before health() calls it
# SHIFTED.  Like the span window, this is a property of the GRATING and not of
# the code: the 24 slices genuinely record different windows, so the threshold
# has to clear the real spread while still catching the ~200-340 A
# misidentifications it exists for.  Measured on setups already gated at 24/24:
#
#     RH2 2023-11-17   max |blue - median|  31.8 A
#     RH2 2024-10-29                        38.4 A
#     RH4 2024-12-28                        34.0 A
#     RH3 2024-10-29 (9 seed slices)        44.6 A   <- only 5 A of headroom
#
# 50.0 stands as the default so nothing that called this before behaves
# differently; pass `shift_tol` for a grating whose real spread runs closer to
# it, and say in that run's notes why.
SHIFT_TOL = 50.0

# NLINE_MIN is a property of the GRATING and the LAMP, not of the code, for the
# same reason the span window is.  20 was set from RM1 and the RH gratings, where
# the FeAr arc against PypeIt's default FeI/ArI/ArII catalogue fits 70-80 lines
# in one ~1450 A slice -- about 5 lines per 100 A.  RM2 sits 1000-3000 A redder,
# where that catalogue runs out: measured on 2023-09-23 B, a setup gated 24/24
# whose 24 dispersions agree to 0.4%, every slice fits only 10-14 lines over
# 1970 A, i.e. 0.6 per 100 A -- eight times thinner.  At 20 that CORRECT setup
# reports zero slices fit to seed a template, and the bootstrap RH3 section 4
# depends on cannot start.
#
# Pass `nline_min` for such a grating and say in that run's notes what the line
# density actually is.  Lowering it does not lower the standard the way it looks:
# a fit through few lines returns a small rms whether or not it is right (the
# "few points, small residual" trap, RH1_PROCEDURE.md), so the check that carries
# the weight at low line counts is dispersion agreement across the 24 slices of
# one setup and between independent nights -- never the rms.

# The 400-900 A window is the physically possible span of ONE slice, and it was
# set from RH1 (~630 A) and RH2 (~764 A).  It is a property of the grating, not
# of the code, so a grating covering more must widen it or every slice is
# rejected as nonsense.  Passed in as `span=(MIN, MAX)`; omitted, the RH1/RH2
# values stand, so nothing that called this before behaves differently.


def _span(span):
    return (SPAN_MIN, SPAN_MAX) if span is None else (float(span[0]), float(span[1]))


def rank(path, quiet=False, span=None, nline_min=None):
    smin, smax = _span(span)
    nlmin = NLINE_MIN if nline_min is None else int(nline_min)
    wv = WaveCalib.from_file(path, chk_version=False)
    good = []
    for i, wf in enumerate(wv.wv_fits):
        spat = int(wv.spat_ids[i])
        if wf is None or wf.pypeitfit is None or wf.wave_soln is None:
            continue
        w = np.asarray(wf.wave_soln).flatten()
        d = np.diff(w)
        wspan = w.max() - w.min()
        nl = 0 if wf.pixel_fit is None else len(wf.pixel_fit)
        med = np.median(d)
        why = []
        if not (wf.rms < RMS_MAX):
            why.append(f"rms {wf.rms:.2f}")
        if nl < nlmin:
            why.append(f"nlines {nl}")
        if not (smin < wspan < smax):
            why.append(f"span {wspan:.0f}A")
        if med <= 0 or (d <= 0).any():
            why.append("non-monotonic")
        elif np.abs(d / med - 1).max() > DISP_TOL:
            why.append(f"disp wobble {100*np.abs(d/med-1).max():.0f}%")
        if not quiet and not why:
            print(f"  PASS slit {i:2d} spat {spat:4d}: rms={wf.rms:.3f} nlines={nl:3d} "
                  f"{w.min():.1f}-{w.max():.1f} A  disp={med:.4f} A/px")
        if not why:
            good.append(dict(idx=i, spat=spat, rms=float(wf.rms), nlines=nl,
                             wmin=float(w.min()), wmax=float(w.max()), disp=float(med)))
    return good


def health(path, span=None, shift_tol=SHIFT_TOL):
    """Did the setup calibrate -- by PypeIt's standards, not the seed criteria.

    Kept separate from rank() on purpose.  rank() is deliberately strict because
    a template seed should be pristine; applying it as a health check calls a
    perfectly good setup a failure.  Measured on 2023-10-16: 4 slices came in at
    rms 0.58-0.69 px against ~0.10 for the rest, which rank() rejects as seeds,
    while PypeIt flagged none of them and reduced all 24.  At 0.289 A/px that is
    a 0.20 A wavelength error -- irrelevant to a throughput continuum.

    So "solved" here means PypeIt produced a fit over a physically possible
    wavelength range and did not flag the slit, which is the number that decides
    whether the setup's science is whole.
    """
    smin, smax = _span(span)
    wv = WaveCalib.from_file(path, chk_version=False)
    n = len(wv.spat_ids)
    solved, rms, blue = [], [], []
    for i, wf in enumerate(wv.wv_fits):
        if wf is None or wf.pypeitfit is None or wf.wave_soln is None:
            continue
        w = np.asarray(wf.wave_soln).flatten()
        if not (smin < w.max() - w.min() < smax) or (np.diff(w) <= 0).any():
            continue
        solved.append(int(wv.spat_ids[i]))
        rms.append(float(wf.rms))
        # Collected HERE, in the same pass and under the same conditions, so it
        # stays index-aligned with `solved`.  It used to be built in a separate
        # comprehension with a weaker condition (any slit with a wave_soln),
        # which made the two lists different lengths -- `zip` then compared each
        # slit against ANOTHER slit's blue end, and the median was taken over
        # unsolved slits.  Measured on RH3 2023-11-12: 4 solved against 20 in the
        # blue list, spat 452 compared against 904.7 A when its own blue end is
        # 8182.5 A, and all four good slices were reported SHIFTED.  Harmless
        # wherever nearly every slit solves (the two lists coincide), which is
        # why it survived RH1, RH2 and RH4; it bites exactly at the holy-grail
        # seed stage, where most slits are expected to fail.
        blue.append(float(w.min()))

    # A slice can keep a normal width, a smooth dispersion and a low RMS while
    # sitting hundreds of A away from its neighbours, because the template was
    # matched to the wrong window.  Width and monotonicity cannot see that; only
    # comparison against the other slices of the same setup can.  Measured on
    # 2023-10-17 C and 2023-10-20 B: 4 slices offset by 201-337 A, which widened
    # their cubes and collapsed the usable range measure_trim could keep.
    if solved:
        medblue = float(np.median(blue))
        shifted = [s for s, b in zip(solved, blue) if abs(b - medblue) > shift_tol]
    else:
        shifted = []

    flagged = []
    slits = glob.glob(os.path.join(os.path.dirname(path), "Slits_*.fits.gz"))
    if slits:
        s = SlitTraceSet.from_file(slits[0], chk_version=False)
        flagged = [int(s.spat_id[i]) for i in np.where(s.mask != 0)[0]]
    return n, solved, flagged, rms, shifted


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    args = sys.argv[1:]
    span = None
    if "--span" in args:
        i = args.index("--span")
        span = (float(args[i + 1]), float(args[i + 2]))
        del args[i:i + 3]
    nline_min = None
    if "--nline-min" in args:
        i = args.index("--nline-min")
        nline_min = int(args[i + 1])
        del args[i:i + 2]
    shift_tol = SHIFT_TOL
    if "--shift-tol" in args:
        i = args.index("--shift-tol")
        shift_tol = float(args[i + 1])
        del args[i:i + 2]
    for p in args:
        print(f"\n=== {p}")
        n, solved, flagged, rms, shifted = health(p, span=span,
                                                  shift_tol=shift_tol)
        ok = [s for s in solved if s not in flagged and s not in shifted]
        print(f"  {len(ok)} usable of {n} slices"
              + (f", {len(flagged)} flagged by PypeIt" if flagged else "")
              + (f", {len(shifted)} SHIFTED vs neighbours {shifted}" if shifted else "")
              + (f"  (rms median {np.median(rms):.3f}, max {max(rms):.3f} px)"
                 if rms else ""))
        g = rank(p, span=span, nline_min=nline_min)
        print(f"  {len(g)} of those meet the stricter seed criteria for "
              f"template building")
        if g:
            best = min(g, key=lambda r: (r["rms"], -r["nlines"]))
            print(f"  best seed = spat {best['spat']} "
                  f"({best['wmin']:.1f}-{best['wmax']:.1f} A, rms {best['rms']:.3f})")
