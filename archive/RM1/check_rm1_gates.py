#!/usr/bin/env python
"""One line per RM1 configuration: did its wavelength calibration actually work.

Reads every WaveCalib under the RM1 run and reports the numbers that decide
whether a setup is fit to reduce.  Written because RM1's broken configurations
do not announce themselves: PypeIt logged no `Not enough useful IDs`, wrote all
ten calibration products including the flat field, and returned 0, on setups
where half the slices carry a dispersion 30% away from the truth.

`disp spread` is the column to read first.  The 24 slices of one setup are the
same grating at the same angle, so their dispersions must agree to a fraction of
a percent; a spread of tens of percent means slices were fitted to lines that
are not the lines they were told they are.

Both windows are RM1's, not the defaults: a slice spans ~1453 A (against the
400-900 A that pick_seed_slits inherits from RH1/RH2), and RM1's genuine
slice-to-slice blue-end spread reaches 84 A from the median (against the 50 A
SHIFTED tolerance set from RH2/RH4).

Usage:  python check_rm1_gates.py [--span 1200 1700] [--shift-tol 100]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pypeit_test"))
from pick_seed_slits import rank, health              # noqa: E402
from gate_rm1_template import shift_tol_for           # noqa: E402
from pypeit.wavecalib import WaveCalib                # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--span", type=float, nargs=2, default=[1200.0, 1700.0])
    ap.add_argument("--shift-tol", type=float, default=None,
                    help="override; otherwise picked per SLICER, because "
                         "RM1's genuine blue-end spread is 79-84 A on "
                         "Large/2x2 and 106 A on Small/1x1")
    args = ap.parse_args()
    span = (args.span[0], args.span[1])

    print(f"{'config':<24}{'cenw':>6}{'slicer':>7}{'tmpl':>32}"
          f"{'ok/n':>8}{'seed':>6}{'rmsMed':>8}{'rmsMax':>8}"
          f"{'disp':>9}{'spread':>8}{'blueUnion':>10}")
    for wc in sorted(glob.glob(os.path.join(HERE, "*", "pypeit_run",
                                            "keck_kcrm_*", "Calibrations",
                                            "WaveCalib_*.fits"))):
        sdir = os.path.dirname(os.path.dirname(wc))
        parts = os.path.relpath(wc, HERE).split(os.sep)
        cfg = f"{parts[0]}/{parts[2].replace('keck_kcrm_','')}"
        pf = glob.glob(os.path.join(sdir, "*.pypeit"))
        txt = open(pf[0]).read() if pf else ""
        mc = re.search(r"cenwave:\s*([0-9.]+)", txt)
        cw = int(round(float(mc.group(1)) / 10.0) * 10) if mc else 0
        md = re.search(r"decker:\s*(\S+)", txt)
        dk = md.group(1) if md else "?"
        mt = re.search(r"reid_arxiv\s*=\s*(\S+)", txt)
        tm = os.path.basename(mt.group(1)) if mt else "(default)"

        try:
            stol = (args.shift_tol if args.shift_tol is not None
                    else shift_tol_for(dk))
            n, solved, flagged, rms, shifted = health(wc, span=span,
                                                      shift_tol=stol)
            ok = [s for s in solved if s not in flagged and s not in shifted]
            good = rank(wc, quiet=True, span=span)
            wv = WaveCalib.from_file(wc, chk_version=False)
            disp, blue = [], []
            for wf in wv.wv_fits:
                if wf is None or wf.wave_soln is None:
                    continue
                w = np.asarray(wf.wave_soln).ravel()
                disp.append(np.median(np.diff(w)))
                blue.append(w.min())
            disp, blue = np.array(disp), np.array(blue)
            dmed = np.median(disp)
            spread = 100 * (disp.max() - disp.min()) / dmed
            union = blue.max() - blue.min()
            flag = "" if spread < 2.0 else "   <-- SLICES DISAGREE"
            print(f"{cfg:<24}{cw:>6}{dk:>7}{tm:>32}"
                  f"{len(ok):>4}/{n:<3}{len(good):>6}"
                  f"{np.median(rms) if len(rms) else np.nan:>8.3f}"
                  f"{np.max(rms) if len(rms) else np.nan:>8.3f}"
                  f"{dmed:>9.4f}{spread:>7.1f}%{union:>9.0f}A{flag}")
        except Exception as e:
            print(f"{cfg:<24}{cw:>6}{dk:>7}{tm:>32}   unreadable: "
                  f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
