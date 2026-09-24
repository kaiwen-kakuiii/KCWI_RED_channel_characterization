#!/usr/bin/env python
"""Flag pixels where the sky model diverged, so they are excluded downstream.

The b-spline sky fit can run away at slit/detector corners, producing sky model
values orders of magnitude above any real sky line while SCIIMG is perfectly
normal. BPMMASK marks these pixels good, so a coadd would use them.

This sets the USER bit (bit 11) in BPMMASK for pixels where

    SKYMODEL  >  factor * (99.99th percentile of SKYMODEL over good pixels)

Real sky lines sit near that percentile; blowups are orders above it.
The edit is reversible -- clear the USER bit to undo.

Usage:
    python mask_sky_blowups.py [--factor 10] [--dry-run] spec2d_*.fits ...
"""
import sys
import numpy as np
from astropy.io import fits
from pypeit.images.imagebitmask import ImageBitMask


def process(path, factor=10.0, dry_run=False):
    bm = ImageBitMask()
    with fits.open(path, mode='readonly' if dry_run else 'update') as h:
        names = {x.name: i for i, x in enumerate(h)}
        isky = [i for n, i in names.items() if n.endswith('SKYMODEL')][0]
        ibpm = [i for n, i in names.items() if n.endswith('BPMMASK')][0]
        sky, bpm = h[isky].data, h[ibpm].data
        good = bpm == 0
        s = sky[good]
        s = s[np.isfinite(s)]
        if s.size == 0:
            print(f"  {path}: no good pixels"); return 0
        thresh = factor * np.percentile(s, 99.99)
        bad = good & np.isfinite(sky) & (sky > thresh)
        n = int(bad.sum())
        if n and not dry_run:
            bpm[bad] = bm.turn_on(bpm[bad], 'USER')
            h[ibpm].data = bpm
            h.flush()
        print(f"  {path.split('/')[-1][:52]:52s} thresh={thresh:10.3g}  flagged={n}")
        return n


if __name__ == '__main__':
    args = [a for a in sys.argv[1:]]
    dry = '--dry-run' in args
    if dry: args.remove('--dry-run')
    factor = 10.0
    if '--factor' in args:
        i = args.index('--factor'); factor = float(args[i+1]); del args[i:i+2]
    if not args: sys.exit(__doc__)
    tot = sum(process(p, factor, dry) for p in args)
    print(f"\ntotal pixels flagged: {tot}" + ("  (dry run, nothing written)" if dry else ""))
