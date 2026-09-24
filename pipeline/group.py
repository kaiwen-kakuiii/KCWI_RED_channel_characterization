#!/usr/bin/env python
"""Split a setup's spec2d files into cubes and write one .coadd3d each.

Splits by star, then by a pointing change larger than 3 arcsec, then by an
airmass spread above the grating file's max_airmass_spread. A cube with fewer
than min_frames is not written: one spec2d with combine=True is rejected by
PypeIt, and combine=False writes a one-degree WCS.

    python pipeline/group.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def load_spec2d(paths, need_rate):
    from astropy.io import fits
    rows = []
    for path in paths:
        header = fits.getheader(path)
        row = {
            "spec2d": os.path.basename(path),
            "path": path,
            "star": common.star_name(str(header.get("TARGNAME") or "")),
            "mjd": float(header.get("MJD") or 0.0),
            "airmass": float(header.get("AIRMASS") or 0.0),
            "ra": float(header.get("RAOFF") or 0.0),
            "dec": float(header.get("DECOFF") or 0.0),
            "exptime": float(header.get("EXPTIME") or 0.0),
        }
        if need_rate:
            row["rate"] = frame_rate(path)
        rows.append(row)
    return rows


def frame_rate(path):
    """Sky-subtracted counts per second over the whole frame."""
    import numpy as np
    from astropy.io import fits
    try:
        with fits.open(path) as hdu:
            names = [hdu_item.name for hdu_item in hdu]
            sci = next(name for name in names if name.endswith("SCIIMG"))
            sky = next(name for name in names if name.endswith("SKYMODEL"))
            data = (np.asarray(hdu[sci].data, float)
                    - np.asarray(hdu[sky].data, float))
            mask_name = next((name for name in names if name.endswith("BPMMASK")), None)
            good = np.isfinite(data)
            if mask_name is not None:
                good &= np.asarray(hdu[mask_name].data) == 0
            exptime = float(hdu[0].header.get("EXPTIME", 1.0)) or 1.0
        return float(np.nansum(data[good]) / exptime)
    except Exception as exc:
        print(f"  could not read a count rate from {os.path.basename(path)} "
              f"({type(exc).__name__}). The frame is kept.")
        return float("nan")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    path = common.require_file(
        common.pypeit_path(args.grating, args.night, args.setup), "science.py")
    meta = common.setup_meta(open(path).read())
    sdir = os.path.dirname(path)
    spec2d = sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
    if not spec2d:
        raise common.StepError(f"no spec2d under {sdir}/Science. Run science.py.")

    rows = load_spec2d(spec2d, cfg["count_rate_floor"] is not None)
    groups, short, faint = common.group_frames(
        rows, float(cfg["max_airmass_spread"]), cfg["count_rate_floor"],
        cfg["min_frames"], cfg["exptime_floor_s"])
    print(f"{args.grating} {args.night} setup {meta['letter']}: {len(spec2d)} spec2d")
    print(f"  airmass spread {cfg['max_airmass_spread']}  "
          f"pointing tolerance {common.POINTING_TOL_ARCSEC} arcsec  "
          f"min frames {cfg['min_frames']}")
    if short:
        print(f"  dropping {len(short)} frame(s) under {cfg['exptime_floor_s']:.0f}s")
    if faint:
        for frame in faint:
            print(f"  dropping {frame['spec2d']}: count rate {frame['rate']:.3g} "
                  f"is under {cfg['count_rate_floor']} of the star's median")
    if not short and any(frame["exptime"] < float(cfg["exptime_floor_s"]) for frame in rows):
        print(f"  keeping frames under {cfg['exptime_floor_s']:.0f}s: dropping them "
              f"would leave fewer than {cfg['min_frames']}")

    by_star = {}
    for group in groups:
        by_star.setdefault(group["star"], []).append(group)
    written = []
    for star, visits in by_star.items():
        kept = [visit for visit in visits if not visit["refused"]]
        index = 0
        for visit in visits:
            airmass = [frame["airmass"] for frame in visit["frames"]]
            if visit["refused"]:
                print(f"  REFUSED {star}: {len(visit['frames'])} frame(s), "
                      f"under {cfg['min_frames']}. No .coadd3d written.")
                continue
            index += 1
            tag = common.cube_tag(star, args.night, meta["letter"], index, len(kept))
            body = common.render_coadd3d(
                tag, star, args.night, meta["letter"], index, len(kept), visit["frames"])
            dest = os.path.join(sdir, f"{tag}.coadd3d")
            with open(dest, "w") as fh:
                fh.write(body)
            written.append(tag)
            print(f"  {tag}: {len(visit['frames'])} frames, airmass "
                  f"{min(airmass):.2f}-{max(airmass):.2f}  -> {os.path.basename(dest)}")
    if not written:
        raise common.StepError("no cube to build. Every group was under the frame minimum.")
    common.next_step("coadd.py", args.grating, args.night, meta["letter"])


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
