#!/usr/bin/env python
"""Fit one sensitivity function per cube, serially.

The polynomial order comes from the grating file. The algorithm is IR and the
extraction is BOX for every grating. Cubes in one setup are fit one after
another: parallel fits have raced on the telluric cache and died on a
half-written file. Do not start a second sensfunc.py while one is running.

    python pipeline/sensfunc.py --grating RH1 --night 2023-11-08 --setup A
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
    trims = sorted(glob.glob(os.path.join(sdir, "*.trim")))
    if args.tag:
        trims = [path for path in trims if os.path.basename(path).startswith(args.tag)]
    if not trims:
        raise common.StepError(f"no .trim file under {sdir}. Run trim.py.")
    order = int(cfg["polyorder"])
    print(f"{args.grating}: algorithm {common.SENS_ALGORITHM}, "
          f"extr {common.SENS_EXTR}, polyorder {order}")
    print(f"  {len(trims)} cube(s), one at a time")
    log = common.log_path(args.grating, args.night)
    for trim in trims:
        tag = os.path.basename(trim)[:-5]
        blue, red = open(trim).read().split()
        sens_file, sens_fits = common.sens_names(tag, order)
        body = common.render_sens(blue, red, order)
        print(f"{tag}: trim {blue}, {red} -> {sens_fits}")
        if args.dry_run:
            print(body.rstrip())
            continue
        with open(os.path.join(sdir, sens_file), "w") as fh:
            fh.write(body)
        spec1d = sorted(glob.glob(os.path.join(sdir, f"spec1d_{tag}*.fits")))
        if not spec1d:
            raise common.StepError(f"no spec1d for {tag}")
        common.run_logged(
            [common.tool("pypeit_sensfunc"), os.path.basename(spec1d[-1]),
             "-s", sens_file, "-o", sens_fits],
            log, cwd=sdir)
        print(f"  wrote {sens_fits}")
    print("Compare this cube with a good night of the same star before adding")
    print("it to the blaze. A low cube is listed under drop: in the grating file.")
    print("Figures: python pipeline/utils/plot_blaze_summary.py")


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
