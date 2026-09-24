#!/usr/bin/env python
"""Stages [3] and [4] of THROUGHPUT_PROCEDURE.md for every RH2 config.

    measure_trim.py --spec2d Science/  -> trim_std_pixs = blue, red
    pypeit_sensfunc spec1d -s <.sens>  -> zeropoint + throughput

Settings are Kaiwen's call, 2026-08-24: **polyorder 15, algorithm IR**.

    algorithm = IR     joint sensitivity + telluric fit; the procedure doc marks
                       this correct above 7000 A, and RH2 runs 6358-7554 A
    extr      = BOX    aperture sum; the cube path has no object-profile model
                       for optimal extraction to use, so OPT is not meaningful
    polyorder = 15     best of 5/9/15 on the RL set

Note polyorder 15 was tuned on RL configs spanning ~3400 A; an RH2 config spans
~800 A, so the fit is ~4x more flexible per Angstrom here.  This script therefore
reports the same two metrics the RL comparison used -- interior scatter and worst
edge error, both in magnitudes -- so overfitting would be visible rather than
assumed absent.  Interior scatter RISING versus a stiffer fit is the signature.

**Strictly sequential.** Parallel sensfuncs race on the shared telluric cache in
~/.cache/pypeit; one job reads a half-written entry and dies with
`FileNotFoundError .../contents`.

Idempotent: a config whose sens file exists is skipped.
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
ROOT = os.path.dirname(BASE)
TRIM = os.path.join(ROOT, "pypeit_test", "measure_trim.py")

POLYORDER = 15
ALGORITHM = "IR"
SUFFIX = "sens_IR_p15cov"


def sh(cmd, log, cwd):
    with open(log, "a") as f:
        f.write(f"\n$ {' '.join(cmd)}\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                              cwd=cwd).returncode


def summarise(path, tag):
    """Zeropoint fit quality and throughput, read only inside the fitted range.

    Outside SENS_ZEROPOINT_FIT_GPM the polynomial extrapolates and can be wrong
    by magnitudes, so every number here is masked to it.
    """
    with fits.open(path) as hdu:
        d = hdu["SENS"].data
        cols = d.columns.names
        w = np.asarray(d["SENS_WAVE"]).ravel()
        zp = np.asarray(d["SENS_ZEROPOINT"]).ravel()
        zf = np.asarray(d["SENS_ZEROPOINT_FIT"]).ravel()
        gpm = np.asarray(d["SENS_ZEROPOINT_FIT_GPM"]).ravel().astype(bool)
        thr = np.asarray(d["THROUGHPUT"]).ravel() if "THROUGHPUT" in cols else None

    ok = gpm & np.isfinite(w) & (w > 0) & np.isfinite(zp) & np.isfinite(zf)
    if not ok.any():
        print(f"{tag}: sens written but no good fit pixels", flush=True)
        return
    lo, hi = w[ok].min(), w[ok].max()
    resid = zp[ok] - zf[ok]
    ww = w[ok]
    interior = (ww > lo + 300) & (ww < hi - 300)
    edge = (ww < lo + 150) | (ww > hi - 150)
    sc = np.std(resid[interior]) if interior.sum() > 5 else np.nan
    ee = np.max(np.abs(resid[edge])) if edge.sum() else np.nan
    print(f"{tag}: fit range {lo:.0f}-{hi:.0f} A  "
          f"interior scatter {sc:.4f} mag  worst edge error {ee:.3f} mag",
          flush=True)

    if thr is not None:
        probes = [p for p in (6500, 6800, 7100, 7400) if lo < p < hi]
        vals = []
        for p in probes:
            i = np.argmin(np.abs(ww - p))
            t = thr[ok][i]
            vals.append(f"{p}A {100*t:.1f}%")
        if vals:
            print(f"{tag}: throughput  " + "   ".join(vals), flush=True)


for sdir in sorted(glob.glob(os.path.join(BASE, "*", "pypeit_run", "keck_kcrm_*"))):
    s1 = sorted(glob.glob(os.path.join(sdir, "spec1d_*.fits")))
    if not s1:
        continue
    s1 = s1[0]
    night = sdir.split(os.sep)[-3]
    cfg = os.path.basename(sdir).replace("keck_kcrm_", "")
    tag = f"{night}/{cfg}"
    log = os.path.join(sdir, "sens.log")
    # Name after the star actually observed -- RH2 is feige110 on five nights
    # and g191b2b on 2023-09-16.
    star = os.path.basename(s1).replace("spec1d_", "").split(f"_{night}")[0]
    out = f"{star}_{night}_{cfg}_{SUFFIX}.fits"
    if os.path.exists(os.path.join(sdir, out)):
        print(f"{tag}: {out} already present, skipped", flush=True)
        continue

    # ---- [3] trim ----------------------------------------------------------
    # --spec2d restricts the fit to wavelengths ALL slits recorded.  Without it
    # the cube's partial-coverage ends are included, where the coadd fabricates
    # flux by resampling across missing slices -- counts look completely normal
    # there, so no counts threshold can catch it.
    r = subprocess.run([os.path.join(ENV, "python"), TRIM,
                        "--spec2d", os.path.join(sdir, "Science"), s1],
                       capture_output=True, text=True)
    open(log, "a").write(r.stdout + r.stderr)
    m = re.search(r"trim_std_pixs = (\d+),\s*(\d+)", r.stdout)
    if not m:
        print(f"{tag}: TRIM FAILED -- no trim_std_pixs in output, see {log}",
              flush=True)
        continue
    blue, red = int(m.group(1)), int(m.group(2))
    cov = re.search(r"full slice-coverage over .*?: ([\d.]+)-([\d.]+) A", r.stdout)
    print(f"{tag}: trim_std_pixs = {blue}, {red}"
          + (f"   (full coverage {cov.group(1)}-{cov.group(2)} A)" if cov else ""),
          flush=True)

    # ---- [4] sensfunc ------------------------------------------------------
    sfile = os.path.join(sdir, f"{star}_{ALGORITHM}_p{POLYORDER}cov.sens")
    open(sfile, "w").write(
        "[sensfunc]\n"
        f"    algorithm = {ALGORITHM}\n"
        "    extr = BOX\n"
        f"    trim_std_pixs = {blue}, {red}\n"
        f"    polyorder = {POLYORDER}\n")

    print(f"{tag}: sensfunc (algorithm {ALGORITHM}, polyorder {POLYORDER}) ...",
          flush=True)
    rc = sh([os.path.join(ENV, "pypeit_sensfunc"), os.path.basename(s1),
             "-s", os.path.basename(sfile), "-o", out], log, sdir)
    if rc or not os.path.exists(os.path.join(sdir, out)):
        print(f"{tag}: SENSFUNC FAILED rc={rc}, see {log}", flush=True)
        continue
    print(f"{tag}: wrote {out}", flush=True)
    try:
        summarise(os.path.join(sdir, out), tag)
    except Exception as e:
        print(f"{tag}: could not summarise ({type(e).__name__}: {e})", flush=True)

print("ALL SENSFUNCS DONE", flush=True)
