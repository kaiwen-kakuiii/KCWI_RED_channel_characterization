#!/usr/bin/env python
"""Cross-check two independent wavelength solutions where they overlap.

A small fit RMS only says the polynomial passes through the lines the fitter
chose.  It cannot see a misidentification: assign every feature the neighbouring
catalogue line and the fit is just as tight, but every wavelength is wrong by one
line spacing.  Two setups at different central wavelengths, solved independently
from different arcs and different slices, overlap in wavelength -- so if both are
right, the *observed* arc features must land at the same wavelengths.  That is
measured here by cross-correlating the two arc spectra on a common wavelength
grid.

What counts as agreement is set by what the template is for.  full_template uses
the template only to decide which catalogue line each feature is; the wavelength
solution it reports is then fit to catalogue wavelengths, not to the template.
So the template's own scale has to be good to a fraction of the *line spacing*
(here ~12 A), not to a fraction of a pixel.  The threshold below is 25% of the
median spacing; anything under that identifies the same lines.

Comparing the two solutions' `wave_fit` arrays instead would prove nothing: those
are catalogue values, identical whenever both picked the same line, whatever
pixel each put it at.

Usage:
    python check_seed_agreement.py <WaveCalib_a.fits> <spat_a> <WaveCalib_b.fits> <spat_b>
"""
import sys

import numpy as np
from pypeit.wavecalib import WaveCalib


def seed(path, spat):
    wv = WaveCalib.from_file(path, chk_version=False)
    i = int(np.where(wv.spat_ids == int(spat))[0][0])
    wf = wv.wv_fits[i]
    return (np.asarray(wf.wave_soln).flatten(),
            np.asarray(wv.arc_spectra[:, i]).astype(float),
            np.asarray(wf.wave_fit).flatten())


def main(pa, sa, pb, sb):
    wa, fa, la = seed(pa, sa)
    wb, fb, lb = seed(pb, sb)

    lo, hi = max(wa.min(), wb.min()), min(wa.max(), wb.max())
    if hi - lo <= 0:
        sys.exit(f"no overlap: {wa.min():.0f}-{wa.max():.0f} vs {wb.min():.0f}-{wb.max():.0f}")

    step = min(np.median(np.diff(wa)), np.median(np.diff(wb)))
    grid = np.arange(lo, hi, step)

    def prep(w, f):
        g = np.interp(grid, w, f)
        # running-median continuum: the two setups sit at different points of the
        # blaze, and that broad shape would otherwise drive the correlation.
        k = max(int(30.0 / step) | 1, 5)
        pad = np.pad(g, k // 2, mode="edge")
        cont = np.array([np.median(pad[i:i + k]) for i in range(g.size)])
        g = g - cont
        return g / (g.std() or 1.0)

    ga, gb = prep(wa, fa), prep(wb, fb)

    nlag = int(8.0 / step)
    lags = np.arange(-nlag, nlag + 1)
    cc = np.array([np.sum(ga * np.roll(gb, k)) for k in lags]) / ga.size
    j = int(np.argmax(cc))
    k = float(lags[j])
    if 0 < j < len(cc) - 1:                       # parabolic peak refinement
        y0, y1, y2 = cc[j - 1], cc[j], cc[j + 1]
        k += 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    shift = k * step

    both = np.sort(np.concatenate([la[(la > lo) & (la < hi)], lb[(lb > lo) & (lb < hi)]]))
    spacing = float(np.median(np.diff(both))) if both.size > 2 else 12.0
    tol = 0.25 * spacing

    print(f"overlap {lo:.1f}-{hi:.1f} A ({hi - lo:.0f} A), "
          f"median catalogue line spacing {spacing:.1f} A")
    print(f"cross-correlation peak {cc.max():.3f} at {shift:+.3f} A "
          f"({shift / step:+.1f} binned px, {abs(shift) / spacing:.2f} line spacings)")
    print(f"zero-shift correlation {cc[np.argmin(np.abs(lags))]:.3f}")

    ok = abs(shift) < tol
    print(f"\n{'AGREE' if ok else 'DISAGREE'}: offset {abs(shift):.2f} A vs "
          f"tolerance {tol:.2f} A (25% of a line spacing).")
    print("  Line identification is unaffected at this level -- full_template refits\n"
          "  to catalogue wavelengths, so the template's own scale need only be\n"
          "  close enough to name each line."
          if ok else
          "  At this level the template could name the wrong catalogue line.")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
