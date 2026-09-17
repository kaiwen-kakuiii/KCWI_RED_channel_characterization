#!/usr/bin/env python
"""Re-run pypeit_sensfunc with measure_trim.py's cut removed, for comparison.

The *_p15cov / *_p5cov products cannot be un-trimmed after the fact: the
zeropoint fit never saw the pixels measure_trim.py removed, so the polynomial
does not exist there.  Recovering that range means re-fitting from the same
spec1d with the trim set to its minimum.

`trim_std_pixs = 1, 1` is that minimum -- pypeit does `trim_gpm[blue:-red]`, so
a 0 gives an empty slice and masks the whole spectrum.  One pixel per end is the
smallest legal cut, and it leaves the reduction's full delivered range in the
fit.  Everything else -- algorithm, extraction, polynomial order -- is held at
the settled value for that grating, so the only difference between a *notrim
product and its *cov counterpart is our trim.

What survives in a *notrim product is then pypeit's own narrowing:
SENS_ZEROPOINT_FIT_GPM, the sigma-rejection inside the zeropoint fit.

These products are for comparison only.  The partial-slice-coverage zones they
restore are the ones measured to be biased by up to ~0.2 mag with normal-looking
counts (THROUGHPUT_PROCEDURE, step 3), which is why the trim exists.

Usage:
    python run_notrim_sens.py [--jobs 3] [--dry-run] [--force]
"""
import argparse
import concurrent.futures as cf
import glob
import os
import subprocess
import sys
import time

import plot_blaze_summary as P

SENSFUNC = "/opt/miniconda3/envs/pypeit/bin/pypeit_sensfunc"

# polyorder is settled per grating and must not change here: the comparison is
# about the trim, so every other knob is held.
POLYORDER = {"IR_p15cov": 15, "IR_p5cov": 5}

HEADER = """# IR sensfunc, polyorder={order}, measure_trim.py's cut REMOVED.
# trim_std_pixs = 1, 1 is the minimum pypeit accepts (it does trim_gpm[b:-r], so
# a 0 masks the whole spectrum).  The reduction's full delivered range goes into
# the fit and only pypeit's own sigma-rejection narrows it.
#
# FOR COMPARISON WITH {cov} ONLY.  This restores the
# partial-slice-coverage zones, which are biased by up to ~0.2 mag with
# normal-looking counts -- that bias is the reason the trim exists.
[sensfunc]
    algorithm = IR
    extr = BOX
    trim_std_pixs = 1, 1
    polyorder = {order}
"""


def jobs():
    """One job per config the blaze analysis actually uses."""
    out = []
    for g in P.ORDER:
        cov = P.SUFFIX[g]
        new = cov.replace("cov", "notrim")
        pat = os.path.join(P.ROOT, f"{g} pypeit run", "*", "pypeit_run",
                           "keck_kcrm_*", f"*_sens_{cov}.fits")
        for sens in sorted(glob.glob(pat)):
            if "order_study" in sens:
                continue
            sdir = os.path.dirname(sens)
            tag = os.path.basename(sens).split("_sens_")[0]
            if tag in P.DROP:                 # non-photometric, never plotted
                continue
            s1 = os.path.join(sdir, f"spec1d_{tag}.fits")
            if not os.path.exists(s1):
                continue
            out.append(dict(grating=g, dir=sdir, tag=tag, cov=cov, new=new,
                            order=POLYORDER[cov],
                            spec1d=os.path.basename(s1),
                            par=f"{tag}_notrim.sens",
                            out=f"{tag}_sens_{new}.fits"))
    return out


def run(j, force=False):
    dest = os.path.join(j["dir"], j["out"])
    if os.path.exists(dest) and not force:
        return j, "skip", 0.0
    with open(os.path.join(j["dir"], j["par"]), "w") as fh:
        fh.write(HEADER.format(order=j["order"], cov=j["cov"]))
    t0 = time.time()
    r = subprocess.run([SENSFUNC, j["spec1d"], "-s", j["par"], "-o", j["out"]],
                       cwd=j["dir"], capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0 or not os.path.exists(dest):
        return j, "FAIL: " + r.stderr.strip().splitlines()[-1][:160], dt
    return j, "ok", dt


def main():
    ap = argparse.ArgumentParser()
    # The telluric cache races when several sensfuncs start cold (one job found
    # a half-written entry and died); 3-wide on a warm cache has been fine, and
    # anything that does fail is retried serially below.
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    js = jobs()
    print(f"{len(js)} configs", flush=True)
    if a.dry_run:
        for j in js:
            print(f"  {j['grating']:4} {j['tag']:34} -> {j['out']}")
        return

    done, failed = 0, []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for j, st, dt in ex.map(lambda x: run(x, a.force), js):
            done += 1
            print(f"[{done:2d}/{len(js)}] {j['grating']:4} {j['tag']:34} "
                  f"{st:8} {dt:5.0f}s", flush=True)
            if st.startswith("FAIL"):
                failed.append(j)

    for j in failed:                          # serial retry, per the cache race
        j2, st, dt = run(j, force=True)
        print(f"[retry] {j['tag']:34} {st} {dt:.0f}s", flush=True)

    bad = [j["tag"] for j in js
           if not os.path.exists(os.path.join(j["dir"], j["out"]))]
    print(f"missing after retry: {bad if bad else 'none'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
