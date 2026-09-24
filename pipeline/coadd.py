#!/usr/bin/env python
"""Build one datacube per .coadd3d file, with combine = True.

After each cube, checks the WCS scale. A degree-scale CDELT is the
single-frame failure and this step stops rather than hand the cube to
extraction.

    python pipeline/coadd.py --grating RH1 --night 2023-11-08 --setup A
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline import common


def check_cube(path, boxcar):
    """Return False when the WCS scale is the broken single-frame kind."""
    import numpy as np
    from astropy.io import fits
    from scipy import ndimage

    with fits.open(path) as hdu:
        header = hdu[1].header if len(hdu) > 1 else hdu[0].header
        cube = np.asarray(hdu[1].data, float)
    cd = [abs(float(header.get(f"CDELT{i}", np.nan))) for i in (1, 2)]
    naxis = [int(header.get(f"NAXIS{i}", 0)) for i in (1, 2)]
    broken = any(value > common.CDELT_MAX_DEG for value in cd if np.isfinite(value))
    spaxel = cd[1] * 3600.0 if np.isfinite(cd[1]) else 1.0
    white = np.nansum(np.where(np.isfinite(cube), cube, 0.0), axis=0)
    total = white.sum()
    local = np.nan
    radius = (boxcar / spaxel) if spaxel > 0 else 2.5
    if total > 0 and radius > 0:
        peak = _star_peak(white, radius, ndimage)
        yy, xx = np.mgrid[0:white.shape[0], 0:white.shape[1]]
        distance = np.hypot(yy - peak[0], xx - peak[1])
        surround = white[distance <= 3 * radius].sum()
        if surround > 0:
            local = white[distance <= radius].sum() / surround
    field = [n * c * 3600.0 for n, c in zip(naxis, cd)]
    print(f"  cube {cube.shape}  CDELT {cd[0]:.3g} x {cd[1]:.3g} deg  "
          f"spaxel {spaxel:.3f} arcsec  field {field[0]:.0f} x {field[1]:.0f} arcsec")
    if np.isfinite(local):
        print(f"  {100 * local:.1f}% of the flux within 3 apertures sits inside "
              f"the {boxcar} arcsec aperture")
    if broken:
        print("  WCS scale is a degree per spaxel. Refusing to extract.")
        return False
    if np.isfinite(local) and local < 0.60:
        print("  WARNING: the star is not compact. Look at the whitelight image "
              "before extracting.")
    return True


def _star_peak(white, radius, ndimage):
    import numpy as np
    box = max(2, int(round(radius / 2.0)))
    smoothed = ndimage.uniform_filter(white, size=box, mode="nearest")
    margin = int(np.ceil(radius))
    ok = np.zeros(white.shape, bool)
    if 2 * margin < min(white.shape):
        ok[margin:-margin, margin:-margin] = True
    else:
        ok[1:-1, 1:-1] = True
    if not ok.any():
        return np.unravel_index(np.argmax(smoothed), smoothed.shape)
    return np.unravel_index(np.argmax(np.where(ok, smoothed, -np.inf)), smoothed.shape)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grating", required=True)
    parser.add_argument("--night", required=True)
    parser.add_argument("--setup", required=True)
    parser.add_argument("--tag", default=None, help="one cube tag; default is every .coadd3d")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = common.load_config(args.grating)
    sdir = common.setup_dir(args.grating, args.night, args.setup)
    files = sorted(glob.glob(os.path.join(sdir, "*.coadd3d")))
    if args.tag:
        files = [path for path in files if os.path.basename(path).startswith(args.tag)]
    if not files:
        raise common.StepError(f"no .coadd3d under {sdir}. Run group.py.")
    print(f"{args.grating} {args.night} setup {common.letter_of(args.setup)}: "
          f"{len(files)} cube(s), combine = True, boxcar {cfg['boxcar_arcsec']} arcsec")
    log = common.log_path(args.grating, args.night)
    for coadd3d in files:
        print(os.path.basename(coadd3d))
        common.run_logged(
            [common.tool("pypeit_coadd_datacube"), os.path.basename(coadd3d), "-o"],
            log, cwd=sdir, dry=args.dry_run)
        if args.dry_run:
            continue
        cube = os.path.join(sdir, os.path.basename(coadd3d).replace(".coadd3d", ".fits"))
        if not os.path.isfile(cube):
            raise common.StepError(f"coadd did not write {cube}")
        if not check_cube(cube, float(cfg["boxcar_arcsec"])):
            raise common.StepError(f"{os.path.basename(cube)} has a broken WCS. Not extracting.")
    common.next_step("extract.py", args.grating, args.night, common.letter_of(args.setup))


if __name__ == "__main__":
    try:
        main()
    except common.StepError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
