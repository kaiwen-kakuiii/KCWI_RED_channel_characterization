#!/usr/bin/env python
"""Audit a wavelength scale against the night sky, at any wavelength.

Same purpose as `check_skyline_wavelengths.py`: the arc cannot audit itself.  If
the template names each feature as the neighbouring catalogue line, every slice
inherits the same offset and every fit is still tight.  Night-sky emission is in
the science exposure, not the calibration, and nothing in the reduction ever fit
to it, so it is an independent ruler.

**Why this exists beside that script.**  `check_skyline_wavelengths.py` centroids
four named lines hardcoded at 6300-6922 A.  RH4 runs 9418-10402 A, where none of
them exist, and where there is no [OI] singlet at all -- everything is OH, in
molecular bands dense enough that individual centroids scatter.  So this script
matches the *whole* OH forest at once by cross-correlation rather than measuring
single lines: with ~70 lines in band, a common-mode shift shows up far more
sharply than any one blend could show it.

**This matters most for RH4.**  RH1 and RH2 could gate a template on a different
night from the one that seeded it.  RH4 is a single night and a single
configuration, so its gate re-solves the arc the template came from -- a
self-consistency check.  The sky is the one ruler in this data set that the
template did not help draw.

PypeIt's OH line lists are in **vacuum**, verified against the air values in
`check_skyline_wavelengths.py`: OH 6498.729 air -> 6500.525 vac, and
`OH_GMOS_lines.dat` carries 6500.521.  So no air/vacuum conversion is applied
here, and applying one would fake a ~2.7 A offset at 10000 A -- larger than the
error being looked for.

Usage:
    python check_skylines_xcorr.py <spec2d_*.fits> [more ...]
    python check_skylines_xcorr.py --lines OH_R24000 --max-shift 15 <spec2d.fits>
"""
import argparse
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

LISTDIR = ("/opt/miniconda3/envs/pypeit/lib/python3.13/site-packages/pypeit/"
           "data/arc_lines/lists")
# Densest OH catalogue covering 8000-25000 A; 107 lines in RH4's 9550-10400.
DEFAULT_LIST = "OH_NIRSPEC_Y"


def sky_spectrum(path, step=0.2):
    """1D sky spectrum binned by WAVELENGTH, not by column.

    Every pixel of an IFU spec2d carries its own wavelength -- the 24 slices each
    see a slightly different window -- so the sky model is accumulated into
    wavelength bins across the whole frame rather than collapsed along a column.
    """
    with fits.open(path) as hdu:
        det = [e.name for e in hdu if e.name.endswith("SKYMODEL")][0].split("-")[0]
        sky = hdu[f"{det}-SKYMODEL"].data
        wave = hdu[f"{det}-WAVEIMG"].data
    good = np.isfinite(wave) & (wave > 100) & np.isfinite(sky)
    lo, hi = float(wave[good].min()), float(wave[good].max())
    edges = np.arange(lo, hi, step)
    idx = np.digitize(wave[good], edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    prof = np.zeros(len(edges) - 1)
    cnt = np.zeros(len(edges) - 1)
    np.add.at(prof, idx[ok], sky[good][ok])
    np.add.at(cnt, idx[ok], 1.0)
    w = 0.5 * (edges[:-1] + edges[1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        prof = np.where(cnt > 0, prof / np.maximum(cnt, 1), np.nan)
    return w, prof, (lo, hi)


def model_spectrum(w, lines, amps, fwhm):
    """Catalogue lines as Gaussians on the same grid, so the two are comparable."""
    sig = fwhm / 2.3548
    m = np.zeros_like(w)
    for lw, a in zip(lines, amps):
        near = np.abs(w - lw) < 5 * sig
        if near.any():
            m[near] += a * np.exp(-0.5 * ((w[near] - lw) / sig) ** 2)
    return m


def detrend(x):
    """Remove the continuum so the correlation is driven by lines, not by shape."""
    x = np.where(np.isfinite(x), x, 0.0)
    n = len(x)
    k = max(11, (n // 40) | 1)
    pad = np.pad(x, k // 2, mode="edge")
    base = np.convolve(pad, np.ones(k) / k, mode="valid")[:n]
    y = x - base
    s = np.std(y)
    return y / s if s > 0 else y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec2d", nargs="+")
    ap.add_argument("--lines", default=DEFAULT_LIST,
                    help="PypeIt OH line list name (vacuum), without _lines.dat")
    ap.add_argument("--fwhm", type=float, default=None,
                    help="line FWHM in A; default is 5 pixels of the observed grid")
    ap.add_argument("--max-shift", type=float, default=10.0,
                    help="largest shift searched, in A")
    ap.add_argument("--step", type=float, default=0.2)
    ap.add_argument("--tol", type=float, default=None,
                    help="offset that counts as agreement, in A.  Default is 25%% "
                         "of the median catalogue line spacing in band, the same "
                         "rule check_seed_agreement.py uses -- the template only "
                         "has to name each line correctly")
    args = ap.parse_args()

    path = os.path.join(LISTDIR, args.lines + "_lines.dat")
    if not os.path.exists(path):
        sys.exit(f"no such line list: {path}")
    t = Table.read(path, format="ascii.fixed_width", comment="#")
    cw = np.asarray(t["wave"], float)
    ca = np.asarray(t["amplitude"], float) if "amplitude" in t.colnames \
        else np.ones(len(cw))

    for p in args.spec2d:
        print(f"\n=== {os.path.basename(p)}")
        w, obs, (lo, hi) = sky_spectrum(p, step=args.step)
        sel = (cw > lo) & (cw < hi)
        if sel.sum() < 5:
            print(f"  only {sel.sum()} catalogue lines in {lo:.0f}-{hi:.0f} A "
                  f"-- cannot audit with {args.lines}")
            continue
        spacing = float(np.median(np.diff(np.sort(cw[sel]))))
        tol = args.tol if args.tol is not None else 0.25 * spacing
        fwhm = args.fwhm if args.fwhm is not None else 5 * args.step * 2.5

        mod = model_spectrum(w, cw[sel], ca[sel], fwhm)
        o, m = detrend(obs), detrend(mod)

        nlag = int(round(args.max_shift / args.step))
        lags = np.arange(-nlag, nlag + 1)
        cc = np.array([np.mean(o * np.roll(m, L)) for L in lags])
        best = int(np.argmax(cc))
        shift = lags[best] * args.step          # A the MODEL moved to match data
        # The absolute correlation is LOW here and that is expected, not a
        # warning: OH catalogue intensities do not predict observed line
        # strengths well, and the sky is faint in a short standard exposure.
        # What matters is whether the peak stands above the correlation's own
        # noise, so it is reported in sigma against the far wings of the lag
        # curve.  Measured on RH4's three stacked feige110 frames: correlation
        # 0.141, but 14.4 sigma.
        nw = max(5, len(cc) // 8)
        wings = np.concatenate([cc[:nw], cc[-nw:]])
        sig = ((cc[best] - wings.mean()) / wings.std()) if wings.std() > 0 else np.nan
        # A positive shift means the catalogue had to move red to match the data,
        # i.e. the reduction's wavelengths are too BLUE by that amount.
        print(f"  sky spans {lo:.1f}-{hi:.1f} A, {int(sel.sum())} {args.lines} "
              f"lines in band, median spacing {spacing:.2f} A")
        print(f"  cross-correlation peak {cc[best]:.3f} at {shift:+.2f} A "
              f"({shift/spacing:+.2f} line spacings); zero-shift {cc[nlag]:.3f}")
        print(f"  peak significance {sig:.1f} sigma above the lag curve's wings"
              + ("   *** peak is not significant -- offset is not measured ***"
                 if np.isfinite(sig) and sig < 4 else ""))
        verdict = "AGREE" if abs(shift) <= tol else "OFFSET"
        print(f"  {verdict}: |{shift:+.2f}| A against a {tol:.2f} A tolerance "
              f"(25% of a line spacing)")
        if abs(shift) > tol:
            print("  The arc cannot see this: a common-mode misidentification "
                  "leaves every per-slice RMS small.  Check the template before "
                  "trusting the wavelengths.")


if __name__ == "__main__":
    main()
