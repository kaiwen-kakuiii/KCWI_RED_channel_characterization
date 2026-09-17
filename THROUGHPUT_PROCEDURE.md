# KCRM RL throughput from feige110 — procedure

How a reduced spec2d becomes a throughput curve, what every setting means, and
what is known to break. Written 2026-08-20.

Scope of the run this documents: 7 nights, 9 configs, 27 science frames of
feige110 through the RL grating. One config (2025-05-25) has no sensfunc — see
*Known problems*.

---

## Vocabulary

Terms used throughout, defined once:

| term | meaning |
|---|---|
| **spec2d** | reduced 2D detector image for one exposure: 24 IFU slices side by side |
| **slice** | one of the 24 strips the image slicer cuts the field into |
| **spaxel** | one spatial pixel of the datacube; carries a full spectrum |
| **pixel** (in a spec1d) | one wavelength bin of the extracted 1D spectrum |
| **whitelight image** | the cube collapsed over all wavelengths — a picture of the field |
| **zeropoint** | magnitude describing instrument sensitivity at a wavelength. Higher = more sensitive. Ours run ~20.5 |
| **throughput** | fraction of photons arriving at the telescope that get recorded. Ours ~29.5% |
| **slice wavelength window** | the wavelength range one slice records; each slice's window is shifted ~few px (total spread ~465 A for RL) |
| **full slice-coverage interval** | wavelengths recorded by ALL 24 slices = largest window start .. smallest window end |
| **partial-coverage zone** | wavelengths near the cube ends recorded by only SOME slices (~465 A at each end) |

---

## The pipeline

```
  spec2d files (one per exposure)
       |
       v
  [1] pypeit_coadd_datacube            ~9 min per config
       |    resample the 24 slices of every exposure onto one 3D grid
       v
  datacube  (30 x 23 spaxels x ~2300 wavelength bins)  + whitelight image
       |
       v
  [2] pypeit_extract_datacube          ~1 min
       |    sum the light in a circle around the star, at every wavelength
       v
  spec1d  (counts vs wavelength)
       |
       v
  [3] pypeit_test/measure_trim.py      instant
       |    decide which wavelengths are usable:
       |    counts criterion AND full slice-coverage
       v
  trim_std_pixs = blue, red
       |
       v
  [4] pypeit_sensfunc                  ~2 min
       |    counts / known feige110 spectrum -> zeropoint -> smooth fit
       v
  zeropoint + THROUGHPUT + QA plots
```

---

## [1] Coadd

Write a `.coadd3d` file listing that config's spec2d files:

```
[rdx]
    spectrograph = keck_kcrm
[reduce]
    [[cube]]
        combine = True
        output_filename = feige110_<night>_<cfg>.fits
        save_whitelight = True

spec2d read
filename
Science/spec2d_....fits
spec2d end
```

Run `pypeit_coadd_datacube feige110.coadd3d -o`.

Each exposure's slices land at slightly different wavelengths and sky positions;
this step puts them on a common grid. `save_whitelight` also writes the collapsed
image, which step 2 needs to find the star.

---

## [2] Extraction

`pypeit_extract_datacube <cube>.fits -o`

Fits a 2D Gaussian to the whitelight image to get the star's centre and width
sigma, then sums flux inside a circle of radius **4 sigma** at every wavelength.

4 sigma encloses 99.97% of a Gaussian, so the aperture adapts to seeing while
always capturing the same *fraction* of the light. Measured radii across these
nights: 2.00"-2.72". This is the right choice for throughput — a fixed aperture
would lose more light on a blurry night and understate the instrument.

Output `spec1d` carries two extractions: `BOX_*` (aperture sum) and `OPT_*`
(optimally weighted). **Use BOX** — there is no object-profile model in the cube
path for optimal extraction to use.

---

## [3] Trim — which wavelengths are usable

This step is not part of pypeit. It exists because of step 4's limitation.

The RL grating's efficiency collapses below ~5550 A: feige110 delivers 2-5 counts
there against ~3900 mid-band. Those are real measurements of near-zero
sensitivity, but the measured zeropoint falls ~1.3 mag in 83 A, and no smooth
polynomial can follow that. They must be excluded by *defining the fit's
wavelength range*, not left for the fitter to reject.

```
pixel   0 .............. 221 ===================== 2057 ......... 2310
wave 5142 A ........... 5550 A ================ 8936 A ....... 9402 A
     |<-- below thresh ->|<--- kept, fitted --->|<- below thresh ->|
       counts 2-50            counts ~3900        falling / spikes
     trim_std_pixs =  221   ,                                  253
```

**Method** (`pypeit_test/measure_trim.py`): keep the **largest contiguous block**
of pixels whose counts exceed 20% of the 95th-percentile counts.

*Contiguous* matters. Taking simply the first and last pixel above threshold lets
an isolated spike anchor the range: 2024-07-06 C has one cosmic ray at 9175 A
reading 3013 counts, separated by an 8-pixel gap (491 counts) from the real
spectrum. Including it dragged the fit boundary 239 A into dead spectrum and was
the sole cause of that config's three symptoms — worst fluxed/model ratio (1.83),
non-monotonic response to polynomial order, and inflated edge errors. Fixing the
trim brought its ratio to 1.13, the best of the set.

The threshold (20%) is not delicate: counts jump 4.8 -> 52 -> 3690 between 5151
and 5800 A, so anything from 5% to 50% lands in nearly the same place.

**The red value must be >= 1.** pypeit does `trim_gpm[blue:-red]`, so 0 yields an
empty slice and masks the entire spectrum.

### Second criterion: full slice-coverage (`--spec2d`)

The counts criterion is not enough. Each slice mirror feeds the grating at a
slightly different angle, so each slice records a slightly different wavelength
window — same exposure, different λ range per slice:

```
slice 12   5987 ────────────────────────────── 9778
slice 10     5997 ──────────────────────────────── 9797
  ...            (each next slice starts a bit redder)
slice 23                  6452 ──────────────────────────────── 10246
                              |========================|
                              full slice-coverage interval
           ^^^^^^^^^^^^^^^^^^^                          ^^^^^^^^^^^^^
           partial-blue zone (465 A)                    partial-red zone
```

(numbers: 2024-12-28 B). The cube spans the union, so near its ends only some
slices contribute. Mechanism (verified on 2024-12-28 B): the star's light is
spread over several slices (there: 59% / 20% / 12% in slits 12 / 13 / 11, whose
windows start 5987 / 6196 / 6209 A — below 6200 A a third of the star has no
detector data). The coadd's resampling then FILLS the cube spaxels of the
missing slices by smearing flux from neighboring covered slices (weighted-mean
resampling is only flux-conserving where coverage is complete), so the border
region is counted roughly twice: the cube holds 1.40x the 6500-7000-relative
flux at 6010-6100 A where the detector slits show only ~1.15x — a fabricated
~22%, matching the observed +0.23 mag zeropoint excess. Ruled out by direct
measurement: the extraction's own coverage renormalization (`flxscl`, +0.6%),
the residual-sky term (<1%), and the relative scale correction (not applied).
The red-edge throughput "decline" (33% -> 21%) was the same artifact biased the
other way. **The counts look completely normal**, and configs agreeing with each
other proves nothing here if they share slice geometry.

`measure_trim.py --spec2d Science/` reads every spec2d frame's WAVEIMG, takes
per-slit min/max good-pixel wavelengths, forms the full-coverage interval
(largest min .. smallest max), intersects over frames (flexure shifts windows
slightly), and restricts the trim to it. Typical cost: ~460 A blue, ~220 A red
beyond what the counts criterion cut — all of it biased data.

Where the counts criterion still matters: the dichroic/grating cutoff can fall
*inside* full coverage (cenwave-7150 configs: cutoff ~5550, full coverage starts
~5558), so neither criterion alone suffices. Final trim = intersection of both,
then a 10-pixel pad inward on each side (guards the cube's resampling-edge
pixels, which also carry normal counts but wrong zeropoint).

---

## [4] Sensfunc

```
pypeit_sensfunc spec1d_<tag>.fits -s feige110_IR_p15cov.sens -o <tag>_sens_IR_p15cov.fits
```

File naming records the settings history: `*_sens_IR.fits` = polyorder 5 counts
trim, `*_p9` / `*_p15` = higher orders counts trim, `*_p15cov` = polyorder 15 +
full slice-coverage trim (current best). Older outputs are kept on disk for
comparison.

with

```
[sensfunc]
    algorithm = IR          joint sensitivity + telluric fit; correct for lambda > 7000 A
    extr = BOX              use the aperture-sum spectrum from step 2
    trim_std_pixs = b, r    from step 3
    polyorder = 15          smoothness of the fit; see below
```

Divides observed counts by the archival feige110 spectrum to get the zeropoint at
each wavelength, then fits a smooth polynomial through those measurements. `IR`
fits a telluric absorption model simultaneously (needs
`TellPCA_3000_26000_R15000.fits`, cached in `~/.cache/pypeit`).

### Choosing polyorder

Tested 5 (pypeit default), 9, 15 on 8 configs. Metrics are the difference between
fitted and measured zeropoint, in magnitudes (0.1 mag ~ 10% flux error):

| | order 5 | order 9 | order 15 |
|---|---|---|---|
| interior scatter (>300 A from either end) | 0.0351 | 0.0319 | 0.0300 |
| worst edge error (within 150 A of an end) | 0.148 | 0.139 | 0.117 |

Order 15 is best on both and shows no sign of fitting noise (that would appear as
interior scatter *rising*). Gains are decelerating, so ~15 is near the useful
limit; orders 20-25 would locate the turnover if wanted.

Caveat: pypeit's rejection is iterative, so a more flexible fit keeps points a
stiffer one discards — the fitted wavelength range shifts slightly with order
(e.g. 5634 -> 5581 -> 5555 A). Configs at different orders do not cover identical
wavelengths.

---

## Reading the output

Use the throughput **only inside the fitted range**, given by
`SENS_ZEROPOINT_FIT_GPM` in the sens file (HDU `SENS`). Outside it the polynomial
extrapolates and can be wrong by magnitudes.

That failure mode is worth recognising: with the fit range extending into the
grating cutoff, the *fitted* zeropoint read 5.5 where the *measured* value was
13.2. Flux calibration divides by sensitivity, so a near-zero denominator sent the
fluxed standard to **813x above** the model at the blue edge, and throughput above
**100%** at the red edge. Both are symptoms of fit-range mismatch, not of bad data.

A second failure mode looks like real signal: **partial-coverage bias** (see
step 3). Sensitivity rising toward a cube edge, or vertical throughput spikes
near the ends of a config's range, mean the cube coadd was fabricating flux by
resampling around missing-slice gaps there. The p15cov products exclude these zones by construction; in
older products the outer ~460 A (blue) / ~220 A (red) of each fit range are
suspect.

Accuracy inside the range: ~3% interior, up to ~12% within 150 A of a boundary.

### Result, order 15 + full slice-coverage trim (p15cov)

```
night / config   slicer   cenwave   fit range (A)   @7000   @7500   @8000   @8500
2023-11-08 B     Medium      8100   6576-9861       27.5%   30.0%   28.6%   29.7%
2023-11-10 B     Medium      7252   5728-9018       28.8%   29.2%   30.5%   29.2%
2024-07-06 B     Medium      7800   6277-9563       26.6%   28.8%   31.0%   29.0%
2024-07-06 C     Medium      7150   5629-8917       27.8%   29.4%   31.1%   27.9%
2024-11-01 C     Large       7150   5578-8912       30.2%   30.3%   28.8%   29.5%
2024-12-04 B     Large       7150   5577-8913       31.1%   28.6%   29.1%   30.2%
2024-12-04 C     Medium      7150   5623-8915       28.6%   31.1%   28.2%   28.1%
2024-12-28 B     Medium      8000   6471-9759       27.4%   29.8%   30.1%   31.4%

  @7000 A   mean 28.5%   std 1.4%
  @7500 A   mean 29.7%   std 0.8%
  @8000 A   mean 29.7%   std 1.1%
  @8500 A   mean 29.4%   std 1.0%
```

Throughput is flat at **~29.5% from 7500-8500 A** with ~1% scatter across 14
months, two slicers and four central wavelengths. Interior values match the
superseded counts-trim run to 0.2% — the coverage trim removes the corrupt ends
without touching good data. Union of fit ranges still covers 5580-9860 A.

### Superseded: order 15, counts trim only (kept for comparison)

Files `*_sens_IR_p15.fits`. The ends of these fit ranges (up to ~460 A blue,
~220 A red) lie in partial-coverage zones and are biased; interior is fine.

```
night / config   fit range (A)   @7000   @7500   @8000   @8500
2023-11-08 B     6111-10089      27.6%   30.2%   28.7%   29.4%
2023-11-10 B     5555-9247       29.2%   29.5%   30.2%   28.8%
2024-07-06 B     5799-9581       26.6%   28.7%   31.0%   29.1%
2024-07-06 C     5550-8931       27.7%   29.4%   30.9%   27.8%
2024-11-01 C     5554-9138       30.0%   30.3%   29.0%   30.0%
2024-12-04 B     5543-8927       31.1%   28.7%   29.1%   30.0%
2024-12-04 C     5549-8929       28.5%   31.2%   28.1%   28.1%
2024-12-28 B     5998-9982       27.3%   29.7%   30.2%   31.3%
```

---

## Known problems

**Single-frame configs.** `combine = True` with one input trips a pypeit type
check (`check_inputs`: "should all be of type 'list'"). Setting `combine = False`
gets past it but writes a cube with a broken WCS — `CDELT = 1 degree` instead of
~0.00019 — which inflates extracted counts by ~1e21 and makes the sensfunc
meaningless. This is why 2025-05-25 (1 frame) has no throughput. Unsolved.

**Parallel sensfuncs race on the telluric cache.** Running 3-wide, one job found
a half-written entry in `~/.cache/pypeit` and died with `FileNotFoundError ...
/contents`. Retry serially; it succeeds.

**2023-11-10 B degraded with polynomial order** (edge error 0.192 -> 0.200 ->
0.268 for orders 5/9/15), opposite to every other config, while counts-trimmed.
Likely explained: its counts-trim range included 154 A of partial-blue and 221 A
of partial-red coverage, and with the coverage trim its order-15 edge error drops
0.268 -> 0.103, in line with the other configs. Order-dependence not retested.

---

## Tools

| file | purpose |
|---|---|
| `pypeit_test/measure_trim.py` | measure `trim_std_pixs` from a spec1d; `--spec2d Science/` adds the full slice-coverage restriction |
| `pypeit_test/plot_sens_edges.py` | one-figure edge diagnostic: counts + zeropoint fit + throughput + edge zooms |
| `pypeit_test/check_skysub.py` | sky-subtraction quality per slice (upstream of this doc) |
| `pypeit_test/find_object_regions.py` | locate the object to set `user_regions`, before any sky fit |
| `pypeit_test/mask_sky_blowups.py` | flag pixels where the sky b-spline diverged |
