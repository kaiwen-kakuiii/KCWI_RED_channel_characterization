# RH3: what differed from the RL, RH1, RH2 and RH4 runs

Only the deltas. The pipeline itself is [PYPEIT.md](PYPEIT.md) and
[THROUGHPUT_PROCEDURE.md](THROUGHPUT_PROCEDURE.md); the missing-template argument
is [RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) and is not repeated here.

Written 2026-08-25, PypeIt 2.0.1, spectrograph `keck_kcrm`.

Scope: **4 nights, 4 usable configs, 9 standard-star frames** — feige110 on three
nights at cenwave 8600 and g191b2b on one night at cenwave 8400, all Medium
slicer, 2x2 binning. Every config solved 24 of 24 slices.

RH3 is the grating where **PypeIt ships a template and it still fails**, where the
line catalogue is adequate to fit but not to identify, and where the largest
error found was not in PypeIt at all but in this project's own polynomial order.
Sections 1, 2 and 6 are the ones worth reading.

---

## 1. A template that exists, and does not fit

`keck_kcwi.py::config_specific_par` maps dispname RH3 to `keck_kcrm_RH3.fits`, and
`keck_kcwi.py:1231` lists RH3 among the supported setups — the only grating in
this project for which that is true. It is still unusable here:

```
shipped keck_kcrm_RH3.fits    8599.4 - 9389.0 A   (2064 px, binspec 2, 0.384 A/px)
these nights, cenwave 8600      ~8136 - 9033 A     ~50% overlap
these nights, cenwave 8400      ~7924 - 8828 A     ~25% overlap
```

`full_template` says so plainly once the failure line is read past:

```
Shift = 271.45; cc = 0.1466      slit 1
Shift = 479.57; cc = 0.2003      slit 2
Shift = 279.15; cc = 0.1639      slit 3
Shift = 898.66; cc = 0.1515      slit 7      <- wanders 600 px between slits
... Not enough useful IDs                       on all 24
```

Adjacent KCRM slices differ by a few pixels, not hundreds. A shift jumping 600 px
slit to slit at cc ~ 0.15 is the cross-correlation returning its own noise.
Result: 1 usable slice of 24 at 8600 and 10 at 8400, every one of them at
rms 0.000 through too few points — **0 usable seeds either way**.

The shipped template also implies 0.384 A/px against the 0.4045 A/px measured
here, which is a second, independent sign it was built from a different setup.

**A supported grating is not the same as a supported configuration.** Nothing in
PypeIt warns that the arxiv template is centred 400 A away; the failure looks
exactly like a bad template.

---

## 2. A catalogue can be enough to FIT and not enough to IDENTIFY

This cost a wasted run and is the most transferable thing in this file.

RH4's rule is to check catalogue coverage before blaming the data. Done here, it
looks reassuring:

```
band 8136-9033 A     FeI 6   ArI 14   ArII 3   = 23 lines, median spacing 25 A
detected, middle slit, on PypeIt's own processed frames:
    FeAr   32 lines >5 sigma, 14 >20 sigma, 7 >100 sigma
    ThAr   78 lines >5 sigma, 35 >20 sigma, 14 >100 sigma
```

23 catalogue lines against 32 detected is *not* the RH4 pathology (6 against
135), and 23 lines genuinely are enough to constrain the order-4 polynomial
`full_template` fits. On that reading holy-grail was run on FeAr with the default
lamps. It solved **0 of 24**, worse than RH1 (1/24), RH2 (7/24) or RH4 (2/24):

```
slit  spat     rms  nlines     wmin      wmax     span     disp
   1   147   0.000       6  -85773.4   58036.1 143809.5  46.2215
   9   747   0.000       5  -25161.5   24713.9  49875.5 -17.4558
  20  1765   0.000       6  -18252.2   16750.5  35002.7 -11.7349
```

Five to nine lines identified per slit, negative wavelengths, 46 A/px.

**Fitting and identifying are different requirements, and a catalogue can satisfy
one and not the other.** `full_template` is handed the answer — the template says
which catalogue line each feature is — and only then fits, so 23 lines suffice.
Holy-grail has no template and must recognise the *pattern* of line spacings by
brute force; 23 lines over 897 A, median spacing 25 A with a 112 A gap, do not
form patterns distinctive enough to lock onto.

**So a line count that clears RH4's bar says nothing about whether holy-grail can
work.** Count lines against the method that will consume them.

---

## 3. The ThAr route, as on RH4

KCWI takes both lamps and types FeAr `arc`, ThAr `tilt`. `seed_rh3.py` retypes the
three ThAr frames `arc,tilt` and comments the FeAr frames out.

`build_rh4_linelist.py --label RH3 --wmin 7900 --wmax 9100 --amp-min 300` writes
`ThArRH3_lines.dat`:

```
184 lines, 7901.8-9097.3 A, median spacing 4.27 A, max gap 37 A
Th I 158,  Th II 16,  Ar I 10
```

Sized against the frame rather than copied from RH4: 78 detected lines over 835 A
is one per 10.2 A, so a 4.27 A catalogue offers ~2.4 candidates per detected
line. The arc FWHM is 2.1 px = 0.85 A, so 4.27 A spacing is well clear of the
density that invites misidentification. **One band, 7900-9100 A, serves both
cenwaves** — deliberately, so the same list calibrates every night.

`build_rh4_linelist.py` gained `--label` and `--note` for this; run with no
arguments it writes byte-identical output to before, so RH4 is untouched.

### Seeds

```
2024-11-06   12 / 24 usable      <- best of any night in this project
2024-10-29    9 / 24
2023-11-12    1 / 24
2024-10-01    0 / 24             <- see section 4
```

Against RH1's 1/24, RH2's 7/24 and RH4's 2/24, with 58-78 lines identified per
slit instead of 5-9.

### The evidence the solution is real

The same physical slit, solved independently on different nights months apart:

```
slit           2024-10-29   2024-11-06   2023-11-12   agreement
spat   66/67      8136.4       8136.4         -         0.2 A
spat 1468         8156.3       8156.2         -         0.1 A
spat 1919/20      8156.0       8156.0         -         0.1 A
spat 1842         8200.2          -         8200.8      0.6 A
```

A single misidentified line would move a solution by ~4 A at this dispersion.
Agreement at 0.1-0.6 A across independent nights rules that out.

---

## 4. The cenwave-8400 night, and a bootstrap one link longer than RH1's

2024-10-01 seeded 0 of 24, and the cause is in how the night was taken:

```
night        ThAr arc exptime   lines >5 sigma   brightest line
2024-10-01        2.3 s              36            19,367 ADU
2024-10-29       20.0 s              78           151,905 ADU
```

Same lamp, grating, slicer and binning; 8x fainter and less than half the lines.
Its CONTBARS and FLATLAMP frames are short too (40 s against 241 s), so this is
how the night was taken, not one bad frame.

Lowering `sigdetect` to 3 recovers 56 detected lines — comparable to what the good
nights show at 10 sigma — and still seeds nothing. Its four fits are all wrong,
and one is instructive:

```
slit 18  spat 1615   rms 0.294   7 lines   8591.7-9444.5 A   0.4256 A/px
```

That looks plausible in isolation — sensible span, sensible dispersion — and
`health()` counted it as solved. But 8591.7-9444.5 A is **redder than the
cenwave-8600 nights**, which is impossible for a bluer grating setting. Span and
dispersion checks cannot catch it; only comparison against what the setup's
central wavelength implies can.

### The bootstrap

The cenwave-8600 template covers 82% of this detector, far more than the ~50% at
which the shipped template failed. Gated against it, the night still fails —
0/24, rms median 1.009 px — **but it produces 7 slices that pass the strict seed
criteria**, spanning 7924.5-8827.8 A, exactly where a cenwave-8400 setup should
sit:

```
slit  spat     rms  nlines       range          disp
   1   146   0.137      26   7976.7-8827.8   0.4143 A/px
   5   451   0.179      28   7973.3-8827.5   0.4145
   6   525   0.040      25   7934.7-8787.7   0.4147
  ... 7 total, dispersions agreeing to 0.3%
```

Those seed `keck_kcrm_RH3_8400.fits` (2175 px, 7924.5-8827.4 A, 0.4149 A/px), and
the night then gates **24/24 at rms median 0.160 px**.

**A template that fails to calibrate a setup can still be good enough to seed one
that does.** RH1 and RH2 bootstrapped from holy-grail; RH3 bootstrapped from
another cenwave's template when holy-grail had nothing to give. Recorded because
the intermediate result — "FAIL 0/24" — reads like a dead end and is not one.

I called this night dead before reading past the FAIL line. It was not.

---

## 5. Three defects in shared code, found by RH3 and fixed

### 5.1 `pick_seed_slits.health()` compared each slit against another slit's wavelength

`solved` was filtered on span and monotonicity; `blue` was built in a separate
comprehension with a weaker condition. `zip(solved, blue)` then paired
mismatched entries, and the median was taken over unsolved slits:

```
len(solved) = 4   len(blue) = 20
   spat 452  <-  blue  904.7 A   (as paired by the code)
   spat 452  ->  blue 8182.5 A   (its own blue end)
```

Harmless wherever nearly every slit solves — which is why RH1, RH2 and RH4 never
hit it — and wrong exactly at the seed stage, where most slits are expected to
fail. `blue` is now collected in the same pass under the same conditions.
Verified: RH2's and RH4's gated WaveCalibs report identical numbers before and
after.

### 5.2 The SHIFTED tolerance is a property of the grating

`health()` flagged slices whose blue end sits >50 A from the median. Measured on
setups already gated at 24/24:

```
RH2 2023-11-17   max |blue - median|  31.8 A
RH2 2024-10-29                        38.4 A
RH4 2024-12-28                        34.0 A
RH3 2024-10-29 (9 seed slices)        44.6 A   <- 5 A of headroom
```

RH3's genuine slice-to-slice spread is the largest measured here. Now a
parameter (`shift_tol`, `--shift-tol`) defaulting to 50.0, so nothing existing
moves. RH3's gates did not in the end need it — recorded because the margin is
thin and the next grating may.

### 5.3 One bad pixel truncated two nights, silently

`measure_trim` keeps the largest **contiguous** block above the counts threshold,
which exists to stop an isolated cosmic ray anchoring the range out into dead
spectrum. RH3 hit the same mechanism in reverse. On 2024-10-29, one pixel at
8661.8 A reading **-184.7 counts** split the spectrum:

```
block 1   px    0-1295   8137.2-8661.3 A   length 1296   <- kept, being larger
block 2   px 1297-2168   8662.2-9014.7 A   length  872   <- discarded
```

The discarded region carries 270-520 counts against a 107-count threshold. It is
good data. 2024-11-06 lost 191 A the same way.

`--bridge N` joins runs separated by at most N pixels before the largest is
chosen. `measure_trim`'s own default stays 0, so no earlier grating moves;
`run_rh3_throughput.py` passes 5.

```
             bridge 0           bridge 5          recovered
2024-10-29   8208.5-8657.3 A    8208.5-8964.5 A     +307 A
2024-11-06   8208.4-8773.9 A    8208.4-8965.0 A     +191 A
2023-11-12   8209.2-8965.2 A    unchanged
2024-10-01   8009.4-8776.4 A    unchanged
```

**Nothing warned.** The pipeline reported `sensfunc OK` and produced a plausible
curve; 2024-10-29 read 35.3% at 8600 A, which was pure edge artefact from a fit
boundary that should not have existed. It surfaced only because one night's
fit range looked short next to the others. That the three 8600 nights now share
a kept range to within 0.6 A is itself the check that the fix is right.

---

## 6. Polynomial order: the largest error in this run was ours

**The first RH3 figure was wrong, and it looked fine.**

Three nights of feige110 through the same grating at the same central wavelength
must produce the same *shape*. Night-to-night transparency scales a curve; it does
not move the peak. At the inherited polyorder 15 they disagreed badly:

```
order  A per DOF        peak positions       shape RMS   interior scatter
   15        47   8720, 8665, 8540 A            2.18%         0.0341
    9        75   8725, 8675, 8570              2.01%         0.0356
    7        94   8700, 8670, 8605              1.18%         0.0391
    5       125   8640, 8675, 8655              0.77%         0.0425
    3       188   8645, 8655, 8645              0.44%         0.0434
```

(shape RMS = scatter between nights after dividing out a per-night grey factor,
measured on 8400-8850 A, i.e. >150 A inside every night's fit edges.)

Peaks 180 A apart at order 15, 35 A apart at order 5. **Shape agreement improves
fivefold while the fit to the data degrades 27%.**

### Why order 15 was too flexible here, and why the RL precedent misleads

RL tuned order 15 across ~3400 A — one degree of freedom per 227 A. RH3's fit
range is ~750 A, where order 15 is one per 47 A, five times more flexible. RH2
noted the same ratio problem and kept order 15 because interior scatter stayed
low. **Low scatter is also what a polynomial following noise produces**, and RH3
shows the sharper test: whether independent nights of the same star agree.

Order 5 is chosen. At 125 A per DOF it is still *more* flexible than the RL
tuning, so it cannot be accused of imposing the shape, and it resolves the
disagreement. Order 3 agrees marginally better but at 188 A/DOF approaches a
forced parabola.

Only the *shapes* were fit artefacts. The **levels** were not — 2024-10-29 reads
~2 points high at every order tested, which is real night-to-night transparency
and is what averaging nights is for.

`p3/p7/p9/p15` outputs are kept under each setup's `order_study/`, out of the
default plot glob. `plot_rh1_throughput.py` gained `--sens-glob` for this.

**This is worth revisiting for RH2 and RH4**, which used order 15 over comparable
~700 A ranges. Not done here — RH3's data cannot settle it for another grating,
and re-deriving those runs was out of scope.

---

## 7. Results

Polyorder 5, aperture 3.4 arcsec, trim = counts AND full slice-coverage with
`--bridge 5`.

```
cube                     cenw  airm   fit range (A)    scat    edge  interior   peak
feige110_2023-11-12_B    8600  1.21   8213.6-8965.1  0.0724  0.178     31.0%   32.1% @ 8640
g191b2b_2024-10-01_B     8400  1.21   8009.6-8776.3  0.0459  0.129     29.7%   30.9% @ 8469
feige110_2024-10-29_B    8600  1.31   8213.8-8964.1  0.0354  0.081     32.7%   34.5% @ 8673
feige110_2024-11-06_B    8600  1.10   8213.6-8964.8  0.0197  0.058     30.9%   32.0% @ 8656
```

Plot: [rh3_throughput_curves.png](../results/rh3_throughput_curves.png).

### The peak tracks the grating angle, as RH2 found

```
cenwave 8400   peak 8469 A    +69 A
cenwave 8600   peak 8640 A    +40 A
cenwave 8600   peak 8673 A    +73 A
cenwave 8600   peak 8656 A    +56 A       mean |peak - cenwave| = 60 A
```

RH2 measured 38 A for the same behaviour against RL's 1343 A. RH3 is a blazed
high-dispersion grating behaving like one. **Do not compare RH3 configs at a
fixed wavelength** — the 8400 and 8600 configs sample different points of the
blaze, and at 8600 A they legitimately differ by 3 points.

### Cube geometry is clean

All four come out (~2230, 30, 23) — a 16" x 20" Medium-slicer field, 0.679"
spaxels. Flux inside the 3.4" aperture: 97-98% on the three feige110 cubes, with
outer rings of -0.1% to 1.0%. **RH3 did not inherit RH4's edge-ring pathology**
(18% of signed flux in the outermost row/column there).

g191b2b reads 118.4% inside the aperture with a -20.3% outer ring. Both are the
signed-total artefact RH2.md documents: a sky-subtracted cube has negative
pixels, so the denominator is not the enclosed flux and the ratio is not bounded
by 100%. Read it as "compact", not as photometry.

---

## 8. Checks

### The night sky, from outside the arc/template/catalogue loop

`check_skyline_wavelengths.py` centroids four lines at 6300-6922 A, all below
RH3's range. `check_skylines_xcorr.py` matches the whole OH forest instead. The
default `OH_NIRSPEC_Y` list is the right one here — 55 lines in band, median
spacing 7.97 A, giving a 1.99 A tolerance.

```
frame                   shift    xcorr peak   significance
2024-11-06 #1          -0.80 A      0.708        33.3 sigma
2024-11-06 #2          -0.80 A      0.691        31.8 sigma
2023-11-12 #1          -0.80 A      0.739        30.5 sigma
2023-11-12 #2          -0.80 A      0.722        29.7 sigma
```

Far stronger than RH4's equivalent (peaks 0.14 at 8-14 sigma), because 55 OH
lines sit in band and the sky is bright here.

**Read the -0.80 A honestly: it is a systematic offset, not noise.** Identical on
four frames from nights a year apart, it is ~2 binned pixels. It passes the
tolerance and is far too small to affect a continuum, but it is not zero. It is
not an air/vacuum error — that would be ~+2.3 A at these wavelengths, in the
other direction. Cause not established.

### Against RL, over the whole band rather than a corner of it

RH4 could compare with RL only over ~100 A where RL had 2 configs. RH3's band
sits inside RL's well-supported region — **8 configs, throughput_std ~0.01** —
across the whole range:

```
         RH3 (mean of 3)     RL      diff    RL n_configs
 8400 A       27.9%        29.3%     -1.4         8
 8500 A       30.8%        29.4%     +1.4         8
 8600 A       32.5%        29.7%     +2.8         8
 8700 A       32.7%        30.4%     +2.3         8
 8800 A       31.0%        31.1%     -0.2         8
```

Agreement to ~1.4 points at the band edges, with RH3 sitting 2-3 points **above**
RL at 8600-8700 A. That excess is where RH3's blaze peaks and RL is flat, so its
sign and location are what a blazed grating predicts rather than a discrepancy.

---

## 9. What this run cannot claim

**The cenwave-8600 template generalises; the 8400 one is not tested.**
`keck_kcrm_RH3_8600.fits` was seeded from 2024-10-29 alone and then solved
2023-11-12 and 2024-11-06 at 24/24 — genuine generalisation, stronger than RH1,
RH2 or RH4 could show. `keck_kcrm_RH3_8400.fits` was seeded from the only
cenwave-8400 night, so its gate re-solves the arc it came from, exactly as RH4's
did.

Both templates now carry `SEEDNITE`/`SEEDSPAT`/`NSEED`/`CENWAVE` headers, and
`gate_rh3_template.py` reads provenance from the file rather than inferring it.
It had inferred it, from which WaveCalibs currently pass `rank()` — and after a
night is gated its WaveCalib passes too, so every night looked like a seed source
and 2024-10-29's gate was reported as "independent" when it was reproduction. The
gate now says "cannot be determined" rather than asserting independence it cannot
support.

**One star carries the 8600 measurement.** All three 8600 nights are feige110.
RH2 got an independent flux-calibration check from g191b2b agreeing with
feige110; here g191b2b sits at a different cenwave, so it checks the pipeline but
not the 8600 flux scale at the same grating angle.

**The blue end below ~8010 A is unmeasured**, and the 8400 night's underexposed
arcs are why the run reaches 7924 A only through a bootstrapped template rather
than a directly seeded one. A cenwave-8400 night with normal arc exposures would
settle both.
