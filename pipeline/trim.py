#!/usr/bin/env python
"""Measure trim_std_pixs for one cube and write <tag>.trim.

The cut is the counts run (20 percent of the 95th percentile, bridge pixels
from the grating file) intersected with wavelengths covered by every slice
of every spec2d, then 10 pixels of pad. The red value is at least 1.

    python pipeline/trim.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common

TRIM = os.path.join(common.ROOT, "pipeline", "utils", "measure_trim.py")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--tag", default=None)
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    sdir = common.setup_dir(args.grating, args.night, args.setup)
    cubes = sorted(path for path in glob.glob(os.path.join(sdir, "*.fits"))
                   if os.path.isfile(path.replace(".fits", ".coadd3d")))
    if args.tag:
        cubes = [path for path in cubes if os.path.basename(path).startswith(args.tag)]
    if not cubes:
        raise common.StepError(f"no cube under {sdir}. Run coadd.py.")
    print(f"{args.grating}: bridge {cfg['bridge']} px, "
          f"frac {common.TRIM_FRAC}, pad {common.TRIM_PAD}, full slice coverage")
    science = os.path.join(sdir, "Science")
    for cube in cubes:
        tag = os.path.basename(cube)[:-5]
        spec1d = sorted(glob.glob(os.path.join(sdir, f"spec1d_{tag}*.fits")))
        if not spec1d:
            raise common.StepError(f"no spec1d for {tag}. Run extract.py.")
        proc = subprocess.run(
            [common.tool("python"), TRIM, "--spec2d", science,
             "--bridge", str(int(cfg["bridge"])),
             "--frac", str(common.TRIM_FRAC),
             "--pad", str(common.TRIM_PAD),
             spec1d[-1]],
            capture_output=True, text=True)
        sys.stdout.write(proc.stdout)
        if proc.stderr:
            sys.stdout.write(proc.stderr)
        match = re.search(r"trim_std_pixs = (\d+), (\d+)", proc.stdout)
        if match is None or proc.returncode:
            raise common.StepError(f"measure_trim.py failed for {tag}")
        blue, red = int(match.group(1)), int(match.group(2))
        if red < 1:
            raise common.StepError(f"{tag}: red trim is {red}")
        dest = os.path.join(sdir, f"{tag}.trim")
        with open(dest, "w") as fh:
            fh.write(f"{blue} {red}\n")
        print(f"  {tag}: trim_std_pixs = {blue}, {red}  -> {os.path.basename(dest)}")
    common.next_step("sensfunc.py", args.grating, args.night, common.letter_of(args.setup))


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
