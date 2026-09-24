#!/usr/bin/env python
"""Sort one night into PypeIt setups.

Reads fits/by_night/<grating>/<night>/ and writes
reductions/<grating>/<night>/pypeit_run/. Prints each setup. The next step
is template.py on a setup that has science frames.

    python pipeline/setup.py --grating RH1 --night 2023-11-08
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
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    raw = os.path.join(common.ROOT, "fits", "by_night", args.grating, args.night)
    if not os.path.isdir(raw):
        raise common.StepError(
            f"no raw directory {raw}. download_std_red.py writes "
            f"fits/by_night/<grating>/<night>/.")
    out = os.path.join(common.night_dir(args.grating, args.night), "pypeit_run")
    print(f"{args.grating} {args.night}")
    print(f"  grating file pipeline/config/{args.grating}.yaml")
    print(f"  polyorder {cfg['polyorder']}  bridge {cfg['bridge']}  "
          f"lamps {cfg['lamps']}  boxcar {cfg['boxcar_arcsec']} arcsec")
    print(f"  raw {raw}")

    existing = sorted(glob.glob(os.path.join(out, "keck_kcrm_*", "*.pypeit")))
    if existing:
        print(f"setup files already present ({len(existing)}). Not running pypeit_setup.")
    else:
        os.makedirs(out, exist_ok=True)
        common.run_logged(
            [common.tool("pypeit_setup"), "-s", common.SPECTROGRAPH,
             "-r", raw, "-c", "all", "-d", out],
            common.log_path(args.grating, args.night), dry=args.dry_run)
        existing = sorted(glob.glob(os.path.join(out, "keck_kcrm_*", "*.pypeit")))
        if args.dry_run:
            print("dry run: no setup files written")
            return
        if not existing:
            raise common.StepError(f"pypeit_setup wrote nothing under {out}")

    print()
    for path in existing:
        text = open(path).read()
        meta = common.setup_meta(text)
        sci = common.active(common.parse_rows(text), "science")
        gaps = common.calibration_gaps(common.parse_rows(text))
        letter = meta["letter"]
        cen = common.round_cenwave(meta["cenwave"])
        print(f"setup {letter}: {meta['dispname']}  {meta['decker']}  "
              f"cenwave {meta['cenwave']} -> {cen}  "
              f"binning {meta.get('binning', '?')}  "
              f"{len(sci)} science")
        if meta["dispname"] != args.grating:
            print(f"  dispname is {meta['dispname']}, not {args.grating}. Skip this setup.")
            continue
        if not sci:
            print("  no science frames. Skip this setup.")
            continue
        if gaps:
            print(f"  missing {', '.join(gaps)}. Skip this setup.")
            continue
        print("  ", end="")
        common.next_step("template.py", args.grating, args.night, letter)


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
