#!/usr/bin/env python
"""Locate the object in a KCWI/KCRM IFU frame WITHOUT running sky subtraction.

Uses only a raw science frame plus the Slits calibration, so it can be run
after pypeit_setup + calibrations but before any science reduction. Reports:

  - how many slices carry significant object flux  (informational)
  - where the object sits along the slice          -> gives the user_regions string

Usage:
    python find_object_regions.py <raw_science.fits> <Calibrations_dir> [--spec keck_kcrm]
"""
import sys
import numpy as np
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import rawimage
from pypeit.slittrace import SlitTraceSet
from pypeit.par import pypeitpar
import glob
import os

# Object width (% of slice) at or above which no sky region can be defined.
# Real detections run 26-42%; see the note where this is used.
DEGENERATE_WIDTH = 85.0


def analyse(rawfile, calibdir, spec_name='keck_kcrm', contam_frac=0.02):
    spec = load_spectrograph(spec_name)
    par = spec.default_pypeit_par()['scienceframe']
    # minimal processing: overscan only. No flat, no bias -- we only need morphology.
    for k in ['use_biasimage', 'use_darkimage', 'use_pixelflat', 'use_illumflat',
              'use_specillum', 'mask_cr', 'subtract_scattlight']:
        if k in par['process'].keys():
            par['process'][k] = False
    ri = rawimage.RawImage(rawfile, spec, 1)
    img = ri.process(par['process'])
    sci = img.image

    slits = SlitTraceSet.from_file(glob.glob(os.path.join(calibdir, 'Slits_*.fits.gz'))[0])
    nspat = sci.shape[1]
    xx = np.arange(nspat)[None, :]

    print(f"\n=== {os.path.basename(rawfile)} ===")
    print(f"image {sci.shape}   nslits {slits.nslits}")

    # per-slice: collapse along spectral axis, subtract the slice's own baseline.
    # Contamination is judged against the NOISE, not against total flux: a slice
    # with 0.5% of the flux can still be many sigma above sky and wreck a joint fit.
    prof_all, flux, excess_snr = [], [], []
    for i in range(slits.nslits):
        lo = int(np.median(slits.left_init[:, i]))
        hi = int(np.median(slits.right_init[:, i]))
        sub = sci[:, lo:hi]
        p = np.median(sub, axis=0)                  # along-slice profile
        base = np.percentile(p, 25)                 # sky level in this slice
        resid = p - base
        prof_all.append(resid)
        flux.append(np.clip(resid, 0, None).sum())
        # sigma of the along-slice profile: Poisson on the sky, reduced by the
        # number of spectral rows the median averaged over.
        nrow = sub.shape[0]
        sig = max(np.sqrt(max(base, 1.0) / nrow), 1e-6) * 1.253   # 1.253 = median vs mean
        excess_snr.append(np.max(resid) / sig)
    flux = np.array(flux)
    excess_snr = np.array(excess_snr)
    tot = flux.sum()
    frac = flux / tot if tot > 0 else flux * 0

    # A slice is contaminated if its peak excess is detectable well above noise.
    SNR_THRESH = 20.0
    ncontam = int((excess_snr > SNR_THRESH).sum())
    print(f"\n{'slice':>5} {'flux frac':>10} {'peak S/N':>10}  {'':4}")
    for i, (fr, sn) in enumerate(zip(frac, excess_snr)):
        bar = '#' * int(50 * fr / max(frac.max(), 1e-9))
        mark = ' <-- contaminated' if sn > SNR_THRESH else ''
        print(f"{i:5d} {fr:10.4f} {sn:10.1f}  {bar}{mark}")

    # along-slice position, from the brightest slice
    b = int(np.argmax(flux))
    p = prof_all[b]
    n = p.size
    c = np.clip(p, 0, None)
    cum = np.cumsum(c) / c.sum()
    lo_pct = 100 * np.searchsorted(cum, 0.02) / n
    hi_pct = 100 * np.searchsorted(cum, 0.98) / n
    peak_pct = 100 * int(np.argmax(p)) / n

    print(f"\nbrightest slice: {b}   object spans {lo_pct:.0f}%-{hi_pct:.0f}% along the slice "
          f"(peak at {peak_pct:.0f}%)")
    print(f"slices carrying >{100*contam_frac:.0f}% of flux: {ncontam} of {slits.nslits}")

    # NOTE: peak S/N per slice is reported for information only. Attempts to turn it
    # into a joint_fit-vs-per-slice rule did not validate against known nights
    # (2023-12-17: 3 slices bad; 2023-11-08: 14) -- use STD_CHIS from an actual
    # reduction to make that choice, via check_skysub.py.

    pad = 8.0
    a = max(0.0, lo_pct - pad)
    z = min(100.0, hi_pct + pad)

    # A band this wide leaves no sky to fit.  Emitting it anyway is worse than
    # emitting nothing: PypeIt's read_userregions parses ':0' as [0,0] and '100:'
    # as [res-1,res], so ':0,100:' defines *no sky pixels at all* -- the reduction
    # still runs, check_skysub.py then reports nan sky-zone residuals, and its
    # "sky ok" verdict is comparing against an empty set.  Measured on RH2
    # 2023-07-16, which turned out to have no star in the frame at all.
    #
    # Real detections are nowhere near this: every RH1 night and every good RH2
    # setup gave an object width of 26-42% (:40,66:, :32,74:, :35,65:).
    if (z - a) >= DEGENERATE_WIDTH:
        print(f"\n  NO DETECTION: object would span {a:.0f}%-{z:.0f}% of the slice "
              f"({z - a:.0f}% wide, degenerate at >={DEGENERATE_WIDTH:.0f}%).")
        print(f"    No sky region can be defined, so no user_regions is suggested.")
        print(f"    A real point source occupies ~26-42% of the slice and sits at "
              f"the SAME along-slice position in every slice it touches.")
        print(f"    Brightest slice {b}: peak S/N {excess_snr[b]:.1f}, "
              f"carries {frac[b]:.3f} of the flux; {ncontam}/{slits.nslits} slices "
              f"above S/N {SNR_THRESH:.0f}.")
        print(f"    Check the frame for a source before reducing it.")
        return dict(ncontam=ncontam, nslits=slits.nslits, lo=a, hi=z,
                    peak=peak_pct, detected=False)

    print(f"\n  suggested:  user_regions = :{a:.0f},{z:.0f}:")
    return dict(ncontam=ncontam, nslits=slits.nslits, lo=a, hi=z,
                peak=peak_pct, detected=True)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    analyse(sys.argv[1], sys.argv[2])
