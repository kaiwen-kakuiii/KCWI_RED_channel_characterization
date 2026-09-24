#!/usr/bin/env python
"""Run calibrations and report how many of the 24 slices solved.

A setup that does not come back 24/24 stops. 25 slices is the dead-column
case: add exclude_regions for this night and setup, run template.py again,
then run this step with --rerun.

    python pipeline/calib_gate.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--rerun", action="store_true",
                        help="delete Calibrations/ and run calibrations again")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    path, text = common.read_pypeit(args.grating, args.night, args.setup, "template.py")
    meta = common.setup_meta(text)
    cen = common.round_cenwave(meta["cenwave"])
    entry = common.template_entry(cfg, cen, meta["decker"])
    reid = common.reid_value(args.grating, entry)
    lamps = common.lamps_value(cfg["lamps"])
    values = common.header_values(text)
    if values["reid_arxiv"] != reid:
        raise common.StepError(
            f"reid_arxiv is {values['reid_arxiv']}, the grating file says {reid}. "
            f"Run template.py.")
    if norm(lamps) != norm(values["lamps"]):
        raise common.StepError(
            f"lamps is {values['lamps']}, the grating file says {lamps}. Run template.py.")
    exclude = common.optional_setup(cfg, "exclude_regions", args.night, meta["letter"])
    if norm(exclude) != norm(values["exclude_regions"]):
        raise common.StepError(
            "exclude_regions in the .pypeit does not match the grating file. "
            "Run template.py.")

    tol = common.shift_tol(cfg, meta["decker"])
    span = cfg["slice_span_A"]
    sdir = os.path.dirname(path)
    calib = os.path.join(sdir, "Calibrations")
    print(f"{args.grating} {args.night} setup {meta['letter']}")
    print(f"  span {span[0]:.0f}-{span[1]:.0f} A   shift_tol {tol:.0f} A   "
          f"nline_min {cfg['nline_min']}")

    if args.rerun and os.path.isdir(calib) and not args.dry_run:
        shutil.rmtree(calib)
        print("removed Calibrations/")
    wavecal = sorted(glob.glob(os.path.join(calib, "WaveCalib_*.fits*")))
    if not wavecal:
        common.run_logged(
            [common.tool("run_pypeit"), os.path.basename(path), "-c"],
            common.log_path(args.grating, args.night), cwd=sdir, dry=args.dry_run)
        if args.dry_run:
            return
        wavecal = sorted(glob.glob(os.path.join(calib, "WaveCalib_*.fits*")))
    else:
        print("WaveCalib already present. Reporting it. Pass --rerun to recompute.")
    if not wavecal:
        raise common.StepError(f"no WaveCalib under {calib}")

    report = subprocess.run(
        [common.tool("python"), os.path.join(common.ROOT, "pipeline", "utils", "pick_seed_slits.py"),
         "--span", str(span[0]), str(span[1]),
         "--shift-tol", str(tol),
         "--nline-min", str(int(cfg["nline_min"])),
         wavecal[0]],
        capture_output=True, text=True)
    sys.stdout.write(report.stdout)
    if report.stderr:
        sys.stdout.write(report.stderr)
    if report.returncode:
        raise common.StepError("pick_seed_slits.py failed")
    match = re.search(r"(\d+) usable of (\d+) slices", report.stdout)
    if match is None:
        raise common.StepError("could not read the slice count from pick_seed_slits.py")
    usable, total = int(match.group(1)), int(match.group(2))
    if total != 24:
        raise common.StepError(
            f"{total} slices, not 24. If this is the dead column, add "
            f"exclude_regions for {args.night}/{meta['letter']} to the grating "
            f"file (the value must end with a comma), run template.py, then "
            f"calib_gate.py --rerun.")
    if usable != 24:
        raise common.StepError(
            f"{usable} of 24 slices usable. Do not run regions.py on this setup. "
            f"A shifted slice or the wrong template looks like a successful fit.")
    print("24/24 slices usable")
    common.next_step("regions.py", args.grating, args.night, meta["letter"])


def norm(value):
    return None if value is None else str(value)


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
