#!/usr/bin/env python
"""Edge-effect diagnostic: 1d counts + zeropoint fit + throughput on one figure.

Three full-range panels sharing the wavelength axis, plus two zoom panels on
the blue and red trim edges. Shows exactly what the trim kept, what the fit
rejected, and where the polynomial extrapolates.

    panel 1  BOX_COUNTS from the spec1d (log scale), with the trim threshold
             (20% of the 95th-percentile count level) and the trim edges
    panel 2  SENS_ZEROPOINT points (kept / fit-rejected) + SENS_ZEROPOINT_FIT,
             + the extended zeropoint array (HDU4) that extrapolates past trim
    panel 3  throughput (HDU5), grey where extrapolated, colored inside trim
    bottom   zooms of the two edges: counts (grey fill) + zeropoint + fit

Usage:
    python plot_sens_edges.py <rundir> [suffix]
        rundir  e.g. .../2024-12-28/pypeit_run/keck_kcrm_B
        suffix  sens tag, default 'p15' (uses *_sens_IR_p15.fits)
"""
import sys, os, glob
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(rundir, tag='p15'):
    spec1d = glob.glob(os.path.join(rundir, 'spec1d_*.fits'))[0]
    pat = f'*_sens_IR_{tag}.fits' if tag else '*_sens_IR.fits'
    sensf = glob.glob(os.path.join(rundir, pat))[0]

    # --- spec1d counts (full, untrimmed) ---
    d = fits.open(spec1d)[1].data
    g = lambda k: np.asarray(d[k][0] if d[k].ndim > 1 else d[k]).ravel()
    w, c, iv = g('BOX_WAVE'), g('BOX_COUNTS'), g('BOX_COUNTS_IVAR')
    finite = (iv > 0) & np.isfinite(c) & (w > 0)
    peak = np.percentile(c[finite], 95)
    thresh = 0.20 * peak

    # --- sens file ---
    sn = fits.open(sensf)
    s = sn[2].data
    sw = s['SENS_WAVE'][0].ravel()
    zp = s['SENS_ZEROPOINT'][0].ravel()
    zfit = s['SENS_ZEROPOINT_FIT'][0].ravel()
    gpm = s['SENS_ZEROPOINT_FIT_GPM'][0].ravel().astype(bool)
    ok = sw > 0
    sw, zp, zfit, gpm = sw[ok], zp[ok], zfit[ok], gpm[ok]
    # trim edges = what the sens file actually kept (SENS_WAVE span), so the
    # plot always reflects the trim that produced THIS file
    wtrim_lo, wtrim_hi = sw.min(), sw.max()
    blue = int(np.searchsorted(w, wtrim_lo))
    red = len(w) - 1 - int(np.searchsorted(w, wtrim_hi, side='right') - 1)
    wex = sn[3].data.ravel()          # extended wave  (past the trim)
    zpex = sn[4].data.ravel()         # extended zeropoint (extrapolated)
    thru = sn[5].data.ravel()
    okx = wex > 0
    wex, zpex, thru = wex[okx], zpex[okx], thru[okx]
    inside = (wex >= wtrim_lo) & (wex <= wtrim_hi)

    name = os.path.basename(sensf).replace('.fits', '')
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1.2, 1, 1.1], hspace=0.32, wspace=0.18)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, :], sharex=ax1)
    axB = fig.add_subplot(gs[3, 0])
    axR = fig.add_subplot(gs[3, 1])

    def trim_lines(ax):
        for x in (wtrim_lo, wtrim_hi):
            ax.axvline(x, color='crimson', ls='--', lw=1.2)
        ax.axvspan(w[finite][0], wtrim_lo, color='crimson', alpha=0.06)
        ax.axvspan(wtrim_hi, w[finite][-1], color='crimson', alpha=0.06)

    # panel 1: counts
    ax1.semilogy(w[finite], np.clip(c[finite], 0.3, None), color='0.3', lw=0.6)
    ax1.axhline(thresh, color='darkorange', lw=1.2,
                label=f'trim threshold = 20% of p95 = {thresh:.0f} cts')
    trim_lines(ax1)
    ax1.set_ylabel('BOX_COUNTS')
    ax1.set_title(f'{name}   trim: keep {wtrim_lo:.1f}-{wtrim_hi:.1f} A '
                  f'(trim_std_pixs = {blue}, {red})')
    ax1.legend(loc='lower center', fontsize=9)

    # panel 2: zeropoint
    ax2.plot(wex, zpex, color='0.75', lw=1.0, label='extended zeropoint (HDU4, extrapolated)')
    ax2.plot(sw[gpm], zp[gpm], '.', color='k', ms=2.5, label='zeropoint data, kept by fit')
    ax2.plot(sw[~gpm], zp[~gpm], 'x', color='crimson', ms=3.5, mew=0.8,
             label='zeropoint data, fit-rejected')
    ax2.plot(sw, zfit, color='forestgreen', lw=1.5, label='polynomial fit')
    trim_lines(ax2)
    lo_y = np.percentile(zp[gpm], 1) - 1.5
    hi_y = np.percentile(zp[gpm], 99) + 1.0
    ax2.set_ylim(lo_y, hi_y)
    ax2.set_ylabel('zeropoint (AB mag)')
    ax2.legend(loc='lower center', fontsize=8, ncol=2)

    # panel 3: throughput
    ax3.plot(wex, 100 * thru, color='0.75', lw=1.0, label='throughput, extrapolated zone')
    ax3.plot(wex[inside], 100 * thru[inside], color='steelblue', lw=1.2,
             label='throughput inside trim')
    trim_lines(ax3)
    ax3.set_ylim(0, min(100, 1.3 * 100 * np.nanmax(thru[inside]) if inside.any() else 100))
    ax3.set_ylabel('throughput (%)')
    ax3.set_xlabel('wavelength (A)')
    ax3.legend(loc='lower center', fontsize=9)

    # zoom panels: edge +- pad
    for ax, edge, side in ((axB, wtrim_lo, 'blue'), (axR, wtrim_hi, 'red')):
        pad = 250.0
        sel = (w > edge - pad) & (w < edge + pad) & finite
        axc = ax.twinx()
        axc.fill_between(w[sel], np.clip(c[sel], 0, None), color='0.85', step='mid')
        axc.set_ylabel('counts', color='0.5')
        axc.tick_params(axis='y', colors='0.5')
        zs = (sw > edge - pad) & (sw < edge + pad)
        ax.plot(sw[zs & gpm], zp[zs & gpm], '.', color='k', ms=3)
        ax.plot(sw[zs & ~gpm], zp[zs & ~gpm], 'x', color='crimson', ms=4, mew=0.9)
        ax.plot(sw[zs], zfit[zs], color='forestgreen', lw=1.5)
        xs = (wex > edge - pad) & (wex < edge + pad)
        ax.plot(wex[xs], zpex[xs], color='0.6', lw=1.0)
        ax.axvline(edge, color='crimson', ls='--', lw=1.2)
        if zs.any():
            zz = zp[zs][np.isfinite(zp[zs])]
            ax.set_ylim(max(zz.min() - 0.5, 5), zz.max() + 0.5)
        ax.set_zorder(axc.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.set_xlim(edge - pad, edge + pad)
        ax.set_xlabel('wavelength (A)')
        ax.set_ylabel('zeropoint (AB mag)')
        ax.set_title(f'{side} edge zoom ({edge:.1f} A)', fontsize=10)

    out = os.path.join(rundir, f'{name}_edges.png')
    fig.savefig(out, dpi=130, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'p15')
