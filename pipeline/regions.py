#!/usr/bin/env python
"""Measure the star's sky window and write it for science.py.

Uses the longest science frame above the exposure floor in each
(star, pointing) group. Pointings are rounded to 0.1 arcsec, which separates
the 1.2 arcsec nod whose sky window actually moves. Windows that then agree
to within 3 percent of a slice share one reduction pass.

user_regions is written to object_regions.txt, not to the grating file.

    python pipeline/regions.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common

FINDER = os.path.join(common.ROOT, "pipeline", "utils", "find_object_regions.py")


def measure(raw, calib, py):
    proc = subprocess.run(
        [py, FINDER, raw, calib], capture_output=True, text=True)
    sys.stdout.write(proc.stdout)
    if proc.stderr:
        sys.stdout.write(proc.stderr)
    if "NO DETECTION" in proc.stdout:
        return None
    import re
    match = re.search(r"user_regions = (\S+)", proc.stdout)
    if match is None:
        raise common.StepError("find_object_regions.py did not report a sky window")
    return match.group(1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    path, text = common.read_pypeit(args.grating, args.night, args.setup, "calib_gate.py")
    meta = common.setup_meta(text)
    frames = common.science_frames(text, cfg)
    sdir = os.path.dirname(path)
    calib = os.path.join(sdir, "Calibrations")
    if not os.path.isdir(calib):
        raise common.StepError(f"no {calib}. Run calib_gate.py first.")
    raw = common.raw_dir(text)

    groups = {}
    for frame in frames:
        key = (frame["star"], round(frame["ra"], 1), round(frame["dec"], 1))
        groups.setdefault(key, []).append(frame)

    print(f"{args.grating} {args.night} setup {meta['letter']}: "
          f"{len(groups)} (star, pointing) group(s)")
    measured = []
    py = common.tool("python")
    for key in sorted(groups):
        chosen = max(groups[key], key=lambda frame: frame["exptime"])
        raw_file = os.path.join(raw, chosen["filename"])
        if not os.path.isfile(raw_file):
            raise common.StepError(f"raw frame not found: {raw_file}")
        print(f"  {key[0]} at ({key[1]}, {key[2]}) from {chosen['filename']} "
              f"({chosen['exptime']:.0f}s)")
        region = measure(raw_file, calib, py)
        if region is None:
            print(f"  no star in {chosen['filename']}. That group is not reduced.")
            continue
        measured.append((key, region, [frame["filename"] for frame in groups[key]]))
        print(f"  user_regions = {region}")

    if not measured:
        raise common.StepError("no sky window. Nothing to reduce.")
    clusters = common.cluster_regions(measured)
    regions_path = os.path.join(sdir, "object_regions.txt")
    common.write_regions(regions_path, clusters)
    print(f"{len(measured)} window(s) -> {len(clusters)} reduction pass(es)")
    print(f"wrote {regions_path}")

    # The .pypeit holds one window. science.py rewrites it once per pass.
    values = common.header_values(text)
    rendered = common.render_header(
        text, reid=values["reid_arxiv"], lamps=values["lamps"],
        exclude=values["exclude_regions"], length_range=values["length_range"],
        user_regions=clusters[0]["region"])
    with open(path, "w") as fh:
        fh.write(rendered)
    common.next_step("science.py", args.grating, args.night, meta["letter"])


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
