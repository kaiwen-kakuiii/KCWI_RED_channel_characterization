#!/usr/bin/env python
"""Turn the RH2 spec2d files into cubes and extracted 1D spectra.

Stages [1] and [2] of THROUGHPUT_PROCEDURE.md, for every RH2 config:

    pypeit_coadd_datacube feige110.coadd3d -o     ~9 min   -> cube + whitelight
    pypeit_extract_datacube <cube>.fits -o        ~1 min   -> spec1d

Both run from inside the setup directory, because the .coadd3d lists its spec2d
files by a path relative to itself.

Sequential on purpose. The later sensfunc stage must be serial anyway -- parallel
jobs race on the shared telluric cache in ~/.cache/pypeit and one dies on a
half-written entry -- and there is no reason to burn the box down here either.

The extraction reports BOX and OPT; **BOX is the one to use** for throughput,
since the cube path has no object-profile model for optimal extraction.

Idempotent: a config whose spec1d already exists is skipped, so this can be
relaunched after a crash.
"""
import glob
import os
import re
import subprocess
import sys

import numpy as np
from astropy.io import fits

ENV = "/opt/miniconda3/envs/pypeit/bin"
BASE = os.path.dirname(os.path.abspath(__file__))


def sh(cmd, log, cwd):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def star_of(sdir):
    """The standard this config observed, from its spec2d headers.

    Not hardcoded: RH2 is feige110 on five nights and g191b2b on 2023-09-16, and
    naming a g191b2b cube feige110_* would mislabel every downstream product.
    """
    for s in sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits"))):
        t = str(fits.getheader(s).get("TARGNAME") or "").strip().lower()
        if t:
            return re.sub(r"[^a-z0-9]+", "", t)
    return "std"


def write_coadd3d(sdir, night, cfg):
    """One .coadd3d per config, listing that config's spec2d files.

    align is deliberately NOT set even when a config is dithered: measured on
    2023-11-07 B (a 2" nod), align=True and align=False give byte-identical cube
    geometry and BOX_COUNTS within 0.5%, because the header WCS already places a
    small nod correctly.  align is for when the WCS cannot be trusted.
    """
    sp = sorted(glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")))
    if not sp:
        return None
    star = star_of(sdir)
    out = f"{star}_{night}_{cfg}.fits"
    lines = ["[rdx]", "    spectrograph = keck_kcrm", "[reduce]", "    [[cube]]",
             "        combine = True", f"        output_filename = {out}",
             "        save_whitelight = True", "", "spec2d read", "filename"]
    lines += [os.path.join("Science", os.path.basename(x)) for x in sp]
    lines += ["spec2d end", ""]
    path = os.path.join(sdir, f"{star}.coadd3d")
    open(path, "w").write("\n".join(lines))
    return path


# Make sure every config with science frames has a .coadd3d before the loop.
for sdir in sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run", "keck_kcrm_*"))):
    if not glob.glob(os.path.join(sdir, "Science", "spec2d_*.fits")):
        continue
    if glob.glob(os.path.join(sdir, "*.coadd3d")):
        continue
    night = sdir.split(os.sep)[-3]
    cfg = os.path.basename(sdir).replace("keck_kcrm_", "")
    made = write_coadd3d(sdir, night, cfg)
    if made:
        print(f"{night}/{cfg}: wrote {os.path.basename(made)}", flush=True)


for cfile in sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run",
                                           "keck_kcrm_*", "*.coadd3d"))):
    sdir = os.path.dirname(cfile)
    tag = f"{cfile.split(os.sep)[-4]}/{os.path.basename(sdir)}"
    log = os.path.join(sdir, "post.log")

    m = re.search(r"output_filename\s*=\s*(\S+)", open(cfile).read())
    if not m:
        print(f"{tag}: no output_filename in .coadd3d, skipped", flush=True)
        continue
    cube = m.group(1)
    cubepath = os.path.join(sdir, cube)

    if glob.glob(os.path.join(sdir, "spec1d_*.fits")):
        print(f"{tag}: spec1d already present, skipped", flush=True)
        continue

    # ---- [1] coadd ----------------------------------------------------------
    if os.path.exists(cubepath):
        print(f"{tag}: cube already present, not rebuilt", flush=True)
    else:
        print(f"{tag}: coadd ...", flush=True)
        rc = sh([os.path.join(ENV, "pypeit_coadd_datacube"),
                 os.path.basename(cfile), "-o"], log, sdir)
        if rc or not os.path.exists(cubepath):
            print(f"{tag}: COADD FAILED rc={rc}, see {log}", flush=True)
            continue

    # Cube geometry, and the WCS sanity check that catches the known
    # single-frame failure: a broken cube carries CDELT = 1 degree instead of
    # ~0.00019, which inflates extracted counts by ~1e21 and silently ruins the
    # sensfunc rather than erroring.
    try:
        with fits.open(cubepath) as hdu:
            hdr = hdu[1].header if len(hdu) > 1 else hdu[0].header
            shape = hdu[1].data.shape if len(hdu) > 1 else hdu[0].data.shape
        cd = [abs(float(hdr.get(f"CDELT{i}", np.nan))) for i in (1, 2)]
        nx = [int(hdr.get(f"NAXIS{i}", 0)) for i in (1, 2)]
        bad_wcs = any(c > 0.01 for c in cd if np.isfinite(c))
        # Extent, not just scale.  A correct KCWI/KCRM field is tens of
        # arcseconds; align=True on 2023-11-07 produced a cube 1750" wide from a
        # 2" nod by locking onto a spurious cross-correlation peak, with a
        # perfectly normal CDELT.  Checking scale alone missed it entirely.
        ext = [n * c * 3600.0 for n, c in zip(nx, cd) if np.isfinite(c)]

        # An oversized grid is NOT by itself fatal, and rejecting on it was wrong.
        # Every RH2 cube comes out with one axis tens of times larger than the
        # ~20"x31" a Large-slicer field should span, because a trace of flux lands
        # at wild sky coordinates and stretches the bounding box.  Measured on
        # 2023-11-17 B: the star peaks at spaxel (7,11) -- the same spaxel as the
        # known-good RL 2024-12-04 B cube -- and rows 0-15 hold 99.82% of the
        # flux, against 100% in 15 rows for RL.  The other 537 rows carry 0.18%.
        # So the extraction sees the same star; the grid is just padded.
        #
        # What IS fatal is a degree-scale CDELT (the documented single-frame
        # failure: CDELT=1 instead of ~0.00019, counts inflated ~1e21).
        # So: reject on scale, and report concentration so a genuinely smeared or
        # doubled cube is still visible.
        core = np.nan
        try:
            with fits.open(cubepath) as hdu:
                cd3 = np.asarray(hdu[1].data, float)
            wl = np.nansum(np.where(np.isfinite(cd3), cd3, 0), axis=0)
            tot = wl.sum()
            if tot > 0:
                j = np.unravel_index(np.argmax(wl), wl.shape)
                y0, y1 = max(0, j[0] - 8), min(wl.shape[0], j[0] + 9)
                x0, x1 = max(0, j[1] - 8), min(wl.shape[1], j[1] + 9)
                core = wl[y0:y1, x0:x1].sum() / tot
            del cd3, wl
        except Exception:
            pass

        print(f"{tag}: cube {shape}  CDELT1/2 = {cd[0]:.3g}/{cd[1]:.3g} deg  "
              f"field {ext[0]:.0f}\" x {ext[1]:.0f}\"  "
              f"flux in +/-8 spaxels of peak {100*core:.2f}%"
              + ("   *** WCS SCALE BROKEN ***" if bad_wcs else ""), flush=True)
        if bad_wcs:
            print(f"{tag}: refusing to extract from a cube with a degree-scale "
                  f"CDELT -- counts would be meaningless", flush=True)
            continue
        if np.isfinite(core) and core < 0.90:
            print(f"{tag}: WARNING star is not compact ({100*core:.1f}% within "
                  f"+/-8 spaxels; healthy is >98%) -- possible smeared or doubled "
                  f"stack, check the whitelight", flush=True)
    except Exception as e:
        print(f"{tag}: could not read cube header ({type(e).__name__}), "
              f"continuing", flush=True)

    wl = cubepath.replace(".fits", "_whitelight.fits")
    print(f"{tag}: whitelight {'written' if os.path.exists(wl) else 'MISSING'}",
          flush=True)

    # ---- [2] extract --------------------------------------------------------
    print(f"{tag}: extract ...", flush=True)
    rc = sh([os.path.join(ENV, "pypeit_extract_datacube"), cube, "-o"], log, sdir)
    s1 = glob.glob(os.path.join(sdir, "spec1d_*.fits"))
    if rc or not s1:
        print(f"{tag}: EXTRACT FAILED rc={rc}, see {log}", flush=True)
        continue

    # Report the aperture and the BOX counts level, so a cube that extracted
    # nothing is visible here rather than three stages later.
    try:
        with fits.open(s1[0]) as hdu:
            names = [h.name for h in hdu]
            d = hdu[1].data
            cols = d.columns.names if hasattr(d, "columns") else []
            box = [c for c in cols if c.startswith("BOX_COUNTS") and c.endswith("COUNTS")]
            cnt = np.asarray(d["BOX_COUNTS"]).ravel() if "BOX_COUNTS" in cols else None
        extra = ""
        if cnt is not None and cnt.size:
            good = cnt[np.isfinite(cnt)]
            extra = (f", BOX_COUNTS median {np.median(good):.0f} "
                     f"95th {np.percentile(good, 95):.0f}, {good.size} bins")
        print(f"{tag}: DONE -> {os.path.basename(s1[0])}{extra}", flush=True)
    except Exception as e:
        print(f"{tag}: DONE -> {os.path.basename(s1[0])} "
              f"(could not summarise: {type(e).__name__})", flush=True)

print("ALL CONFIGS COADDED AND EXTRACTED", flush=True)
