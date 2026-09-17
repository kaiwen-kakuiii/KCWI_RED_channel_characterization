#!/usr/bin/env python
"""Measure trim_std_pixs for a spec1d extracted from an IFU datacube.

trim_std_pixs tells pypeit_sensfunc how many pixels to drop from the blue and
red ends of the standard-star spectrum before fitting the sensitivity function.

WHY IT MATTERS
--------------
The ends must be cut back to where the grating actually delivers signal, not
merely to where the counts are non-zero. Inside the RL grating's blue cutoff
(below ~5800 A) feige110 gives 2-5 counts against a peak of ~3900. Those points
carry a real but enormous sensitivity deficit (zeropoint 13 vs 20.5 in the
interior), so pypeit rejects them from the polynomial fit -- and the polynomial
then EXTRAPOLATES into the gap, diverging by several magnitudes. That shows up
as the fluxed standard sitting ~800x above the model at the blue edge, and as
throughput values above 100%.

CRITERION
---------
Keep the contiguous run of pixels whose counts exceed `frac` of the 95th
(`--bridge N` first joins runs separated by <= N pixels, so a single bad
column does not truncate the spectrum at its own position)
percentile (default 0.20), and which have positive inverse variance. Trim
everything outside it, then pad `pad` extra pixels (default 10) inward on each
side. The pad exists because the first/last pixels of the datacube carry
resampling/edge artifacts with NORMAL counts but wrong zeropoint (2024-12-28 B:
first ~15 px sit +0.3 mag above trend at ~6000 counts) -- a counts threshold
cannot see them. The sensfunc fit's own sigma-rejection usually catches them,
but excluding them up front removes that reliance.

FULL SLICE-COVERAGE OPTION (--spec2d)
-------------------------------------
Each IFU slice records a slightly different wavelength window (slice mirrors
feed the grating at different angles; total spread ~465 A for KCRM RL). Near
the cube's wavelength ends only SOME slices contribute ("partial coverage").
The coadd fills the missing slices' spaxels by smearing neighboring slices'
flux (weighted-mean resampling only conserves flux where coverage is
complete), fabricating up to ~20% (2024-12-28 B: cube holds 1.40x the
6500-7000-relative flux at 6010-6100 A vs ~1.15x on the detector = +0.23 mag
vs the 7-config consensus; the extraction's own coverage renormalization
contributes only +0.6%). A counts threshold cannot detect this: the counts
look normal.

With --spec2d <file|dir|glob>, the trim is additionally restricted to the
FULL slice-coverage interval: wavelengths recorded by ALL slits in ALL given
frames (per slit take min/max good-pixel wavelength; full interval = largest
min .. smallest max; intersected over frames, since flexure shifts windows
slightly). Final trim = counts criterion AND full coverage, then pad.

Note the returned values are never 0: pypeit does trim_gpm[blue:-red], so a 0
gives an empty slice and masks the entire spectrum.

Usage:
    python measure_trim.py [--frac 0.2] [--pad 10] [--spec2d Science/] spec1d_*.fits
"""
import sys
import os
import glob
import numpy as np
from astropy.io import fits


def full_coverage(spec2d_arg, verbose=True):
    """FULL slice-coverage interval (lo, hi) in Angstrom over the given frames."""
    if os.path.isdir(spec2d_arg):
        paths = sorted(glob.glob(os.path.join(spec2d_arg, 'spec2d_*.fits')))
    else:
        paths = sorted(glob.glob(spec2d_arg))
    if not paths:
        raise FileNotFoundError(f'no spec2d found for {spec2d_arg}')
    lo, hi = -np.inf, np.inf
    for p in paths:
        h = fits.open(p, memmap=True)
        pick = lambda tag: next(x.data for x in h if x.name.endswith(tag))
        wave, bpm, slits = pick('WAVEIMG'), pick('BPMMASK'), pick('SLITS')
        xx = np.arange(wave.shape[1])[None, :]
        mins, maxs = [], []
        for i in range(len(slits)):
            l = int(np.median(slits['left_init'][i]))
            r = int(np.median(slits['right_init'][i]))
            m = (xx >= l) & (xx <= r) & (wave > 0) & (bpm == 0)
            wv = wave[m]
            mins.append(wv.min())
            maxs.append(wv.max())
        lo, hi = max(lo, max(mins)), min(hi, min(maxs))
    if verbose:
        print(f"full slice-coverage over {len(paths)} frame(s): {lo:.1f}-{hi:.1f} A")
    return lo, hi


def measure(path, frac=0.20, pad=10, coverage=None, verbose=True,
            bridge=0):
    d = fits.open(path)[1].data
    g = lambda k: np.asarray(d[k][0] if d[k].ndim > 1 else d[k]).ravel()
    w, c, iv = g('BOX_WAVE'), g('BOX_COUNTS'), g('BOX_COUNTS_IVAR')
    n = len(c)
    finite = (iv > 0) & np.isfinite(c) & (w > 0)
    if not finite.any():
        print(f"{path}: NO usable pixels"); return None
    peak = np.percentile(c[finite], 95)
    strong = finite & (c > frac * peak)
    if not strong.any():
        print(f"{path}: no pixel above {frac:.0%} of peak"); return None
    idx = np.where(strong)[0]
    # Keep only the LARGEST CONTIGUOUS block of above-threshold pixels. Isolated
    # spikes beyond a gap (cosmic rays, cube artefacts) would otherwise anchor the
    # kept range far out into the cutoff: e.g. 2024-07-06 C has a single pixel at
    # 9175 A reading 3013 counts, separated by an 8-pixel gap (491 counts) from the
    # real spectrum, which dragged the red boundary 130 px past where signal ends.
    splits = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    # ...but a gap of a FEW pixels is a dropout inside the spectrum, not its end.
    # `bridge` joins blocks separated by at most that many pixels before the
    # largest is chosen.  Measured on RH3 2024-10-29: ONE pixel at 8661.8 A
    # reading -184.7 counts split an otherwise continuous spectrum into
    #     px    0-1295  (8137.2-8661.3 A)   and   px 1297-2168  (8662.2-9014.7 A)
    # and the blue block won on length by 424 px, discarding 353 A that carries
    # 270-520 counts against a 107-count threshold.  A single negative outlier is
    # a cosmic ray or bad column; it is not where the spectrum stops.
    #
    # Default 0 keeps the previous behaviour exactly, so RL, RH1, RH2 and RH4
    # results are unchanged unless deliberately recomputed with --bridge.
    if bridge > 0 and len(splits) > 1:
        merged, cur = [], list(splits[0])
        for nxt in splits[1:]:
            if nxt[0] - cur[-1] - 1 <= bridge:
                cur.extend(range(cur[-1] + 1, nxt[0]))   # fill the gap
                cur.extend(nxt)
            else:
                merged.append(np.array(cur)); cur = list(nxt)
        merged.append(np.array(cur))
        splits = merged
    run = max(splits, key=len)
    lo_i, hi_i = run[0], run[-1]
    if coverage is not None:
        inside = np.where((w >= coverage[0]) & (w <= coverage[1]))[0]
        lo_i, hi_i = max(lo_i, inside[0]), min(hi_i, inside[-1])
    # pad inward, but never so much that the kept range collapses
    pad = min(pad, max((hi_i - lo_i - 100) // 2, 0))
    lo, hi = lo_i + pad, hi_i - pad
    blue, red = max(int(lo), 1), max(int(n - 1 - hi), 1)
    if verbose:
        print(f"{path.split('/')[-1][:44]:44s} n={n:5d}  keep {w[lo]:7.1f}-{w[hi]:7.1f} A"
              f"  (dropped {w[0]:6.1f}-{w[lo]:6.1f} blue, {w[hi]:7.1f}-{w[-1]:7.1f} red)"
              f"   trim_std_pixs = {blue}, {red}")
    return blue, red


if __name__ == '__main__':
    a = sys.argv[1:]
    frac = 0.20
    pad = 10
    cov = None
    if '--frac' in a:
        i = a.index('--frac'); frac = float(a[i+1]); del a[i:i+2]
    if '--pad' in a:
        i = a.index('--pad'); pad = int(a[i+1]); del a[i:i+2]
    if '--spec2d' in a:
        i = a.index('--spec2d'); cov = full_coverage(a[i+1]); del a[i:i+2]
    bridge = 0
    if '--bridge' in a:
        i = a.index('--bridge'); bridge = int(a[i+1]); del a[i:i+2]
    for p in a:
        measure(p, frac, pad, coverage=cov, bridge=bridge)
