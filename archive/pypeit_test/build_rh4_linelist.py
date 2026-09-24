#!/usr/bin/env python
"""Build the arc line list RH4 needs, because PypeIt's KCWI default is empty there.

RH4 sits at central wavelength **9950 A**, and that is far enough red that the
line list `keck_kcwi.py` assigns -- `['FeI', 'ArI', 'ArII']` -- has essentially
nothing in it:

    catalogue lines in RH4's 9419-10402  FeI  0    ArI  4    ArII  2   (6 total)
    catalogue lines in RH2's 6400-7600   FeI 25    ArI 40    ArII 13   (78)

`FeI_lines.dat` stops at 9002.0 A and `FeAr_lines.dat` at 8921.9 A.  Both ArII
lines sit in the bluest 100 A, and the two ArI lines covering the red half carry
NIST intensities of 300 and 200 against 35000 for the strong ones.  Six lines
cannot constrain the order-4 polynomial PypeIt fits, so the default configuration
cannot calibrate RH4 no matter what template it is given: `full_template` uses
the template only to decide *which catalogue line* each feature is, and then fits
to catalogue wavelengths (autoid.py, `iterative_fitting(..., line_lists, ...)`).

**The lamp is not the problem -- the catalogue is.**  Measured on the raw frames,
the FeAr lamp puts 114 lines over 100 sigma on the RH4 detector, as many as it
does on RH2.  They are real Fe I and Ar I lines; PypeIt's curated lists simply
stop before 9500 A.

The night carries **ThAr frames as well as FeAr** (KCWI takes both: FeAr is typed
`arc`, ThAr `tilt`), and `ThAr_lines.dat` has 854 lines in 9200-10800 A, 811 NIST-flagged.
So RH4 calibrates off the ThAr frames instead.  See RH4.md.

Two reasons this script exists rather than `lamps = ThAr`:

1. **`lamps = ThAr` alone takes a code path we cannot run.**  `HolyGrail.__init__`
   branches `if 'ThAr' in self._lamps and len(self._lamps) == 1` into
   `run_kdtree()`, which needs a precomputed `ThAr_patterns_poly*_search*.kdtree`.
   That file is not shipped with PypeIt 2.0.1 and the line-list data path has no
   remote host configured (`dataPaths.linelist.host is None`), so it cannot be
   fetched:

       PypeItError: Remote host type None is not supported for package data caching.

   Any other lamp name falls through to `run_brute()`, which needs only the list.

2. **The full list is too dense and too wide.**  17099 lines over 3000-11000 A
   would have brute-force pattern matching search the whole optical for a
   spectrum we already know is centred at 9950 A, and 854 lines over 1600 A is
   one per 1.9 A against an arc line FWHM of about 2 A -- denser than the data
   can resolve, which invites misidentification rather than preventing it.

The output is written into the repo, not into the conda environment, for the
reason `RH1_WAVELENGTH.md` gives for the template: `get_file_path` returns any
path that already resolves, so an absolute path in `lamps` works and an
environment rebuild cannot silently lose the file.

Usage:
    python build_rh4_linelist.py                    # 9200-10800 A, amplitude >= 100
    python build_rh4_linelist.py --amp-min 200      # sparser, if holy-grail misfires
"""
import argparse
import os

import numpy as np
from astropy.table import Table

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = ("/opt/miniconda3/envs/pypeit/lib/python3.13/site-packages/pypeit/"
       "data/arc_lines/lists/ThAr_lines.dat")
# The name PypeIt will rebuild: lamps = <this path without the suffix>
OUT = os.path.join(HERE, "ThArRH4_lines.dat")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=SRC, help="PypeIt's full ThAr line list")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--wmin", type=float, default=9200.0)
    ap.add_argument("--wmax", type=float, default=10800.0,
                    help="band to keep.  Generous either side of RH4's ~850 A "
                         "exposure at cenwave 9950, so the search cannot run out "
                         "of catalogue before it runs out of detector")
    ap.add_argument("--amp-min", type=float, default=100.0,
                    help="NIST relative intensity floor.  100 leaves ~1 line per "
                         "5 A, against ~1 detected line per 6 A measured on the "
                         "raw ThAr frame -- the regime PypeIt expects, where it "
                         "detects a subset of the catalogue rather than the "
                         "catalogue offering more lines than exist")
    ap.add_argument("--label", default="RH4",
                    help="grating this list is for; stamped into the header so a "
                         "list built for another setup is not filed as RH4's")
    ap.add_argument("--note", default=None,
                    help="one line for the header saying why the default "
                         "catalogue was insufficient for that grating")
    ap.add_argument("--keep-unknown", action="store_true",
                    help="keep UNKNWN rows.  Off by default: PypeIt routes those "
                         "through load_unknown_list(), which matches on lamp NAME "
                         "and returns nothing for a name it does not know")
    args = ap.parse_args()

    t = Table.read(args.src, format="ascii.fixed_width", comment="#")
    w = np.asarray(t["wave"], float)
    a = np.asarray(t["amplitude"], float)
    ion = np.asarray(t["ion"])

    keep = (w >= args.wmin) & (w <= args.wmax) & (a >= args.amp_min)
    if not args.keep_unknown:
        keep &= ion != "UNKNWN"
    sub = t[keep]
    sub.sort("wave")
    if len(sub) < 30:
        raise SystemExit(f"only {len(sub)} lines survive -- too few to calibrate "
                         f"with; lower --amp-min")

    with open(args.out, "w") as f:
        f.write(f"# {args.label} arc line list -- ThAr, "
                f"{args.wmin:.0f}-{args.wmax:.0f} A, "
                f"amplitude >= {args.amp_min:.0f}\n")
        f.write(f"# Built by build_rh4_linelist.py from {os.path.basename(args.src)}\n")
        note = args.note or ("PypeIt's default KCWI list (FeI/ArI/ArII) holds 4 "
                             "lines above 9550 A; RH4 is centred at 9950.")
        f.write(f"# {note}\n")
        sub.write(f, format="ascii.fixed_width")

    ws = np.asarray(sub["wave"], float)
    print(f"wrote {args.out}")
    print(f"   {len(sub)} lines, {ws.min():.1f}-{ws.max():.1f} A, "
          f"median spacing {np.median(np.diff(ws)):.2f} A")
    import collections
    print("   by ion:", dict(collections.Counter(np.asarray(sub["ion"])).most_common()))
    print()
    print("   use it with, in the .pypeit file:")
    print(f"       lamps = {args.out[:-len('_lines.dat')]}")


if __name__ == "__main__":
    main()
