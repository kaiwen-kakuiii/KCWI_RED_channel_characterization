#!/usr/bin/env python
"""Write this setup's template, and its line list, into the .pypeit file.

Looks up (central wavelength, slicer) in the grating file. A setting that is
not listed stops here. For a ThAr grating, FeAr arc rows are commented out
and the ThAr rows become arc,tilt. If the template, the lamps, or
exclude_regions changed, Calibrations/ is removed so the old solution cannot
be reused.

    python pipeline/template.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def lamp_of(raw, filename):
    path = os.path.join(raw, filename)
    if not os.path.isfile(path):
        return None
    from astropy.io import fits
    header = fits.getheader(path)
    if str(header.get("IMTYPE")) != "ARCLAMP":
        return None
    if str(header.get("LMP0STAT")) == "1":
        return "FeAr"
    if str(header.get("LMP1STAT")) == "1":
        return "ThAr"
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--dry-run", action="store_true",
                        help="print the parameter block, write nothing")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    path, text = common.read_pypeit(args.grating, args.night, args.setup, "setup.py")
    meta = common.setup_meta(text)
    if meta["dispname"] != args.grating:
        raise common.StepError(
            f"setup {meta['letter']} is {meta['dispname']}, not {args.grating}")
    common.science_frames(text, cfg)

    cen = common.round_cenwave(meta["cenwave"])
    entry = common.template_entry(cfg, cen, meta["decker"])
    reid = common.reid_value(args.grating, entry)
    lamps = common.lamps_value(cfg["lamps"])
    exclude = common.optional_setup(cfg, "exclude_regions", args.night, meta["letter"])
    length = common.optional_setup(cfg, "length_range", args.night, meta["letter"])
    previous = common.header_values(text)

    print(f"{args.grating} {args.night} setup {meta['letter']}")
    print(f"  cenwave {meta['cenwave']} -> {cen} {meta['decker']}")
    print(f"  template {entry}")
    print(f"  reid_arxiv {reid}")
    print(f"  lamps {lamps or '(PypeIt default FeI/ArI/ArII)'}")
    if exclude:
        print(f"  exclude_regions {exclude}")
    if length is not None:
        print(f"  length_range {length}")

    note = None
    if lamps:
        raw = common.raw_dir(text)
        text, note = common.retype_thar(text, lambda name: lamp_of(raw, name))
        print(f"  {note}")

    # A sky window already measured is kept. template.py does not invent one.
    rendered = common.render_header(
        text, reid=reid, lamps=lamps, exclude=exclude, length_range=length,
        user_regions=previous["user_regions"])
    changed = _changed(previous, reid, lamps, exclude, length) or (
        note is not None and note != "already retyped")
    if args.dry_run:
        print(rendered.split("# Setup")[0].rstrip())
        if changed:
            print("dry run: would write the .pypeit and remove Calibrations/")
        return

    with open(path, "w") as fh:
        fh.write(rendered)
    print(f"wrote {path}")
    calib = os.path.join(os.path.dirname(path), "Calibrations")
    if changed and os.path.isdir(calib):
        shutil.rmtree(calib)
        print("removed Calibrations/. The previous solution was keyed by frame "
              "and would have been reused.")
    elif not changed:
        print("wavelength and slit settings already matched the grating file")
    common.next_step("calib_gate.py", args.grating, args.night, meta["letter"])


def _changed(previous, reid, lamps, exclude, length):
    def norm(value):
        return None if value is None else str(value)
    return (norm(previous["reid_arxiv"]) != norm(reid)
            or norm(previous["lamps"]) != norm(lamps)
            or norm(previous["exclude_regions"]) != norm(exclude)
            or norm(previous["length_range"]) != norm(length))


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
