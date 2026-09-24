#!/usr/bin/env python
"""Extract a BOX spectrum at the grating file's boxcar radius, in arcsec.

The radius is in arcsec on purpose. The default aperture is 4 sigma measured
in spaxels, and a spaxel is 1.358 arcsec on Large/2x2 but 0.339 arcsec on
Small/1x1.

    python pipeline/extract.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--tag", default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    sdir = common.setup_dir(args.grating, args.night, args.setup)
    cubes = sorted(glob.glob(os.path.join(sdir, "*.fits")))
    # Sensfunc products and spec1d files also end in .fits. The cube is the
    # file a .coadd3d names.
    cubes = [path for path in cubes
             if os.path.isfile(path.replace(".fits", ".coadd3d"))]
    if args.tag:
        cubes = [path for path in cubes if os.path.basename(path).startswith(args.tag)]
    if not cubes:
        raise common.StepError(f"no cube with a .coadd3d under {sdir}. Run coadd.py.")
    radius = cfg["boxcar_arcsec"]
    print(f"{args.grating}: BOX extraction, boxcar {radius} arcsec")
    log = common.log_path(args.grating, args.night)
    for cube in cubes:
        print(os.path.basename(cube))
        common.run_logged(
            [common.tool("pypeit_extract_datacube"), os.path.basename(cube),
             "-o", "-b", str(radius)],
            log, cwd=sdir, dry=args.dry_run)
    common.next_step("trim.py", args.grating, args.night, common.letter_of(args.setup))


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
