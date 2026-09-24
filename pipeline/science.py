#!/usr/bin/env python
"""Reduce the science frames, one sky window at a time.

Reads object_regions.txt. When two stars or two pointings need different
windows, each pass comments the other science rows out so both still use the
calibrations already built. Frames the exposure floor dropped stay commented.

    python pipeline/science.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def has_spec2d(sdir, filename):
    stem = filename[:-5] if filename.endswith(".fits") else filename
    return bool(glob.glob(os.path.join(sdir, "Science", f"spec2d_{stem}*.fits")))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--rerun", action="store_true",
                        help="run even if spec2d files already exist")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    path, text = common.read_pypeit(args.grating, args.night, args.setup, "regions.py")
    meta = common.setup_meta(text)
    sdir = os.path.dirname(path)
    regions_path = common.require_file(
        os.path.join(sdir, "object_regions.txt"), "regions.py")
    clusters = common.read_regions(regions_path)
    values = common.header_values(text)
    if not values["user_regions"]:
        raise common.StepError("the .pypeit has no user_regions. Run regions.py.")

    print(f"{args.grating} {args.night} setup {meta['letter']}: "
          f"{len(clusters)} pass(es), exposure floor {float(cfg['exptime_floor_s']):.0f}s")
    log = common.log_path(args.grating, args.night)
    for index, cluster in enumerate(clusters, start=1):
        pending = [name for name in cluster["frames"] if not has_spec2d(sdir, name)]
        print(f"pass {index}: user_regions = {cluster['region']}  "
              f"{len(cluster['frames'])} frame(s)")
        if not pending and not args.rerun:
            print("  spec2d already present. Skipping. Pass --rerun to reduce again.")
            continue
        rendered = common.render_header(
            text, reid=values["reid_arxiv"], lamps=values["lamps"],
            exclude=values["exclude_regions"], length_range=values["length_range"],
            user_regions=cluster["region"])
        rendered = common.activate(rendered, set(cluster["frames"]))
        if args.dry_run:
            print(f"  dry run: would reduce {', '.join(cluster['frames'])}")
            continue
        with open(path, "w") as fh:
            fh.write(rendered)
        common.run_logged(
            [common.tool("run_pypeit"), os.path.basename(path)],
            log, cwd=sdir, dry=False)

    if not args.dry_run:
        # Leave every frame that has a sky window active. Floor-dropped frames
        # were never listed in object_regions.txt, so they stay commented.
        kept = [name for cluster in clusters for name in cluster["frames"]]
        final = common.activate(open(path).read(), set(kept))
        with open(path, "w") as fh:
            fh.write(final)
    n2d = len(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
    print(f"{n2d} spec2d file(s) in Science/")
    common.next_step("group.py", args.grating, args.night, meta["letter"])


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
