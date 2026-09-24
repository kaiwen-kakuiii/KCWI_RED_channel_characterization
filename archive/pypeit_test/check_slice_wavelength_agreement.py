#!/usr/bin/env python
"""Find slices whose wavelength solution disagrees with the rest of their setup.

The 24 slices of one setup do not share a wavelength window -- the slicer feeds
the grating at slightly different angles, spreading the windows over ~40 A -- but
they cannot disagree by hundreds.  A slice whose solution is *shifted* keeps a
normal width and a plausible dispersion, so span and monotonicity checks pass it,
and only the comparison against its neighbours exposes it.

Measured on 2023-10-17 C: two slices sat 331 and 337 A blue of the other 22,
with RMS 1.6 and 3.3 px against a median of 0.12.  They widened the cube from
630 A to 946 A, and because measure_trim keeps only the interval covered by ALL
slices, they cut the usable range from ~555 A to 262 A.

Usage:
    python check_slice_wavelength_agreement.py [--tol 50] [<run dir>]
"""
import argparse
import glob
import os
import sys

import numpy as np
from pypeit.wavecalib import WaveCalib


def check(path, tol):
    wv = WaveCalib.from_file(path, chk_version=False)
    rows = []
    for i, wf in enumerate(wv.wv_fits):
        if wf is None or wf.wave_soln is None:
            continue
        w = np.asarray(wf.wave_soln).flatten()
        rows.append((int(wv.spat_ids[i]), float(w.min()), float(w.max()),
                     float(wf.rms)))
    if not rows:
        return None, []
    med = float(np.median([r[1] for r in rows]))
    bad = [r for r in rows if abs(r[1] - med) > tol]
    return (med, len(rows), np.median([r[3] for r in rows])), bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=50.0,
                    help="A a slice's blue end may differ from the setup median")
    ap.add_argument("runs", nargs="?", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "RH1 pypeit run"))
    args = ap.parse_args()

    nbad = 0
    for wc in sorted(glob.glob(os.path.join(args.runs, "*", "pypeit_run",
                                            "keck_kcrm_*", "Calibrations",
                                            "WaveCalib_*.fits"))):
        parts = os.path.relpath(wc, args.runs).split(os.sep)
        tag = f"{parts[0]}/{parts[2]}"
        info, bad = check(wc, args.tol)
        if info is None:
            print(f"{tag:28s} no solved slice")
            continue
        med, n, rms = info
        status = "OK" if not bad else f"{len(bad)} SHIFTED"
        print(f"{tag:28s} {n:2d} slices, blue end {med:7.1f} A, "
              f"rms median {rms:.3f}   {status}")
        for spat, lo, hi, r in bad:
            print(f"      spat {spat:4d}: {lo:7.1f}-{hi:7.1f}  "
                  f"offset {lo - med:+7.1f} A  rms {r:.3f}")
            nbad += 1

    print(f"\n{nbad} shifted slice(s) across all setups"
          + ("" if nbad else " -- every slice agrees with its neighbours"))
    return 1 if nbad else 0


if __name__ == "__main__":
    sys.exit(main())
