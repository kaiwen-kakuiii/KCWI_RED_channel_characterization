# RH4: what differed from the RL, RH1 and RH2 runs

Only the deltas. The pipeline itself is [PYPEIT.md](PYPEIT.md) and
[THROUGHPUT_PROCEDURE.md](THROUGHPUT_PROCEDURE.md); the missing-template argument
is [RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) and is not repeated here. This file is
the problems RH4 raised that those do not cover, and what was done about each.

Written 2026-08-25, PypeIt 2.0.1, spectrograph `keck_kcrm`.

Scope: **1 night, 1 usable config, 5 standard-star frames** — feige110 (3) and
feige34 (2), both on 2024-12-28, Large slicer, 2x2, cenwave 9950. The config
solved 24 of 24 slices.

**RH4 is the first grating where the problem was not the template.** RH1 and RH2
were missing a *reference arc*. RH4 is missing the **line catalogue** — the list
of true wavelengths PypeIt fits to — and no template can substitute for that.
Section 1 is the whole story; everything after it is smaller.

---

## 1. PypeIt's KCWI line list is empty at 9950 A

`keck_kcwi.py::config_specific_par` assigns `lamps = ['FeI', 'ArI', 'ArII']` for
every KCWI/KCRM setup. Counted directly out of the shipped `.dat` files:

| list | file covers | lines in RH4's 9419-10402 A | lines in RH2's 6400-7600 A |
|---|---|---|---|
| `FeI_lines.dat` | 3021.5 - **9002.0** | **0** | 25 |
| `ArI_lines.dat` | 3950.1 - 10883.9 | **4** | 40 |
| `ArII_lines.dat` | 3244.6 - **9511.0** | **2** | 13 |
| `FeAr_lines.dat` | 3021.5 - **8921.9** | **0** | — |
| | | **6 total** | **78 total** |

Six lines, and their distribution is worse than the count:

```
9438.8  ArII      9511.0  ArII      9660.4  ArI      9787.2  ArI
10054.8 ArI (intensity 300)         10335.6 ArI (intensity 200)
|<-- 4 lines in the first 370 A -->|<-- 2 weak lines in the next 610 A -->|
```

Both ArII lines sit in the bluest 100 A, and the two lines carrying the red half
of the range have NIST intensities of 300 and 200 against 35000 for the strong
ones. PypeIt fits an order-4 polynomial with iterative rejection. There is
nothing there to fit.

**This is not fixable with a template.** `full_template` uses the template only to
decide *which catalogue line* each observed feature is; the solution it reports is
then fit to catalogue wavelengths — `autoid.py` ends
`iterative_fitting(obs_spec_i, dets, gd_det, IDs[gd_det], line_lists, ...)`. Build
the best template in the world and it still fits to a catalogue with six lines in
range. This is why the RH1/RH2 recipe cannot simply be re-run for RH4.

### The lamp is fine — it is the catalogue that stops

Worth measuring before concluding the data are bad. Peaks found on a raw arc
frame, same extraction for each, in the middle of the detector:

```
                       >5 sigma   >20 sigma   >100 sigma   brightest line
RH4  FeAr  16.5 s         135        117         114         34893 ADU
RH4  ThAr  16.5 s         149        131         122         34696 ADU
RH2  FeAr  21.5 s         135        130         125         13894 ADU
RH1  FeAr   120  s         75         34           9           949 ADU
```

RH4's arc frames are the *brightest* of the three. Over a hundred real lines are
on the detector; PypeIt's curated lists simply stop before 9500 A.

**Rule, and it is the RH1/RH2 lesson in a new place: when a calibration fails at
an unusual wavelength, check the catalogue's coverage before blaming the data.**

### Fix: calibrate off the ThAr frames

KCWI takes **both** arc lamps. `check_frame_type` types FeAr as `arc` and ThAr as
`tilt`, with the comment "PypeIt is only setup to wavelength calibrate using the
FeAr lamp" — a choice that is right everywhere bluer and wrong here. The night
carries three of each, and `ThAr_lines.dat` holds **854 lines in 9200-10800 A, 811
of them NIST-flagged**.

So `seed_rh4.py` retypes the three ThAr frames `arc,tilt`, comments out the three
FeAr frames, and points `lamps` at a band-limited ThAr list.

### Two reasons the list is rebuilt rather than `lamps = ThAr`

**1. `lamps = ThAr` takes a code path this install cannot run.**
`HolyGrail.__init__` branches

```python
if 'ThAr' in self._lamps and len(self._lamps) == 1:
    self.run_kdtree()
```

which needs a precomputed `ThAr_patterns_poly*_search*.kdtree`. That file is not
shipped with PypeIt 2.0.1, and the line-list data path has no remote host
configured (`dataPaths.linelist.host is None`), so it cannot be fetched:

```
PypeItError: Remote host type None is not supported for package data caching.
```

Any other lamp *name* falls through to `run_brute()`, which needs only the list.

**2. The full list is too wide and too dense.** 17099 lines over 3000-11000 A
would have brute-force pattern matching search the whole optical for a spectrum
whose central wavelength is in the header, and 854 lines over 1600 A is one per
1.9 A against an arc FWHM of about 2 A — denser than the data resolve, which
invites misidentification rather than preventing it.

`build_rh4_linelist.py` therefore writes `ThArRH4_lines.dat`: ThAr, 9200-10800 A,
NIST intensity >= 100, UNKNWN rows dropped.

```
319 lines, 9202.2-10789.3 A, median spacing 3.86 A
Th I 274,  Th II 23,  Ar II 13,  Ar I 9
```

One line per 3.9 A against 1 detected line per 6.6 A on the raw ThAr frame (149
peaks over 5 sigma across 983 A) — the regime PypeIt expects,
where it detects a subset of the catalogue rather than the catalogue offering more
lines than exist.

The UNKNWN rows are dropped because PypeIt routes them through
`load_unknown_list()`, which matches on lamp **name** and returns nothing for a
name it does not recognise.

### The lamp name is an absolute path

`lamps` is not validated against a fixed list (`options['lamps'] = None` in
`WavelengthSolutionPar`), and `load_line_lists` builds `f'{lamp}_lines.dat'` then
calls `get_file_path`, which returns any path that already resolves. So

```ini
lamps = /<repo>/pypeit_test/ThArRH4
```

works, and the list lives in the repo rather than in the conda environment — the
same reasoning [RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) gives for the template,
where an environment rebuild would otherwise lose it silently.

**The template and the line list are one unit.** A run pointed at
`keck_kcrm_RH4.fits` but left on the default lamps fails exactly as if the
template were bad. `run_rh4.py` refuses to start unless both are set.

---

## 2. RH4's slices are wider than the seed filter allowed

`pick_seed_slits.py` rejects any slice whose solution spans outside **400-900 A**,
a window set from RH1 (~630 A) and RH2 (~764 A) to throw out holy-grail's
nonsense fits.

RH4 covers **983-994 A per exposure** at 0.476 A/px. Every slice was rejected, and
the script reported `0 usable of 24` on a WaveCalib that in fact held two good
solutions. Ten minutes were lost to reading that as a failed calibration.

Span is a property of the grating, so it is now a parameter rather than a
constant: `rank()` and `health()` take `span=(MIN, MAX)`, and the CLI takes
`--span 800 1200`. Omitted, the RH1/RH2 values stand, so nothing that called it
before behaves differently — checked against RH2 2023-11-17 B, still 24/24.

**A rejection window derived from one grating is not a property of the code.**

---

## 3. Seeding, and how RH4's gate is weaker than RH2's

Holy-grail against the ThArRH4 list, on the night's only config:

| | usable slices |
|---|---|
| RH1 2023-09-16 B | 1 / 24 |
| RH2 2023-11-17 B | 7 / 24 |
| **RH4 2024-12-28 B** | **2 / 24** |

The run ended non-zero in `flatfield.py illum_profile_spectral`, as on RH1 and
RH2 — WaveCalib is written before the flat field, which is all the seed stage
needs.

The two survivors, and the three "fits" the filter caught:

```
slit  6 (spat  507): rms 0.320  63 lines   9418.6-10401.8 A   0.4761 A/px  <- seed
slit 14 (spat 1302): rms 0.499  61 lines   9426.1-10419.6 A   0.4790 A/px  <- seed
slit  4 (spat  359): rms 0.481  19 lines  10355.4-13591.4 A   1.3633 A/px  <- nonsense
slit  8 (spat  654): rms 0.000   7 lines  10763.4-25484.7 A  -2.1085 A/px  <- nonsense
slit 11 (spat  875): rms 0.524  23 lines   9793.2-10051.7 A   0.1220 A/px  <- nonsense
```

Slit 8 is the RH1 pathology exactly: rms 0.000 through 7 points spanning 15000 A.

The two seeds agree to **0.01 A** over their 976 A overlap (cross-correlation peak
1.001, tolerance 0.77 A). **But read that number for what it is**: RH1 and RH2
compared solutions from *different nights*. These are two slices of the same
exposure, so this is a within-frame consistency check, not independent
confirmation. It is reported here so nobody later mistakes it for the stronger
test.

Because the spans are within 1% of each other (983 vs 994 A) while the RMS differs
0.320 vs 0.499 — the latter sitting on `pick_seed_slits`' 0.5 rejection boundary —
`build_rh4_template.py` takes the **lower-RMS** seed. RH2's "coverage beats RMS"
rule bought 24 A there; 10 A in 985 does not buy anything, so `--span-tol` (2%)
now marks spans that count as equal and the tie goes to RMS.

### The template

```
keck_kcrm_RH4.fits   2064 px, 9418.6-10401.8 A, binspec 2, 0.4761 A/px, monotonic
```

One seed, no stitch — **the good case**. RH1 needed a chain across six central
wavelengths and RH2 needed two chains, one per slicer. RH4's single config gives a
template one exposure wide, which is the shape PypeIt ships: no stitch joint, and
no alignment ambiguity of the kind that cost RH1 four slices
(RH1_PROCEDURE.md problem 3).

### The gate

```
run_pypeit -c returned 0
slices solved, unflagged and unshifted: 24 / 24
rms  median 0.290  max 0.419 binned px
24 meet the stricter seed criteria, spanning 9415.3-10459.9 A
```

2 of 24 to **24 of 24**, and the flat field — which killed the holy-grail run —
completed. Median RMS 0.290 px is 0.138 A at this dispersion, between RH2's Large
(0.058 px) and Medium (0.269 px) templates, and the RL work already established
that a 0.20 A wavelength error is irrelevant to a throughput continuum.

**What this gate cannot claim.** RH2's gate was stronger than RH1's because it ran
the template against a *different night* from the one its seeds came from, which
demonstrates generalisation rather than reproduction. RH4 has one night and one
configuration, so this gate re-solves the arc the template was seeded from. It
still catches what it exists to catch — the unsolved slices that crash the flat
field — and the improvement is real, because identification is what changed while
the fit is to catalogue wavelengths either way. But it is a self-consistency
check. `gate_rh4_template.py` prints that caveat rather than a bare PASS.

The compensating independent check is the night sky, section 5.

---

## 4. Two standards in one configuration

`pypeit_setup` groups on (dispname, decker, binning, cenwave), and **the target is
not part of that key**. RH4's setup B therefore contains both stars:

```
feige110   3 frames   50, 120, 240 s   airmass 1.20-1.22   17:00 UT
feige34    2 frames   180, 180 s       airmass 1.15        04:00 UT (+11 h)
```

**Both existing drivers would have coadded them into a single cube**, and neither
would have said anything:

* `post_rh2.py` reads `TARGNAME` from the *first* spec2d and names the whole config
  after it, then lists every spec2d in one `.coadd3d`.
* `run_rh1_throughput.py` splits on airmass and pointing, but `RAOFF`/`DECOFF` are
  dither offsets and are 0.0 for both stars, and the two sit at airmass 1.20-1.22
  and 1.15 — a spread of 0.07, inside the 0.10 default.

The result would be one cube holding two different stars 11 hours apart,
flux-calibrated against whichever one `pypeit_sensfunc` matched on coordinates.
`run_rh4_throughput.py` splits on `TARGNAME` **first**, before airmass is
considered at all.

### The sky regions differ too, so the reduction splits as well

`find_object_regions.py`, run per star on that star's longest exposure:

```
feige110 (240 s)   user_regions = :33,68:
feige34  (180 s)   user_regions = :38,69:
```

5% of a slice apart — the stars were not acquired at the same place. Reducing both
under one setting would hand 5% of the slice carrying feige110's wings to the sky
model as if it were blank sky, which is the over-subtraction mechanism the RL
skysub work already documented.

So `run_rh4.py` measures the regions per star and, when they differ by more than
`--region-tol` (default 3%), reduces each star in its own pass with the other
star's rows commented out. Calibrations are built once and reused — PypeIt keys
them by frame, not by which science rows are active.

**Two stars is a gain, not a nuisance.** They are an independent check on the flux
calibration, the same cross-check RH2 got from g191b2b agreeing with feige110.

---

## 5. Checks

### The night sky, because the arc cannot audit itself

RH4 needs this more than RH1 or RH2 did, for the reason in section 3: its gate is
self-consistent rather than independent. A common-mode error — every feature named
as the neighbouring catalogue line — leaves every per-slice RMS small and is
invisible to the arc, which defined the scale.

`check_skyline_wavelengths.py` cannot help: it centroids four lines hardcoded at
6300-6922 A, all below RH4's range, and there is no `[OI]` singlet up here at all.
So `check_skylines_xcorr.py` matches the *whole* OH forest at once by
cross-correlation — with ~100 catalogue lines in band, a common-mode shift shows
up far more sharply than any single blend could show it.

PypeIt's OH lists are in **vacuum**, verified rather than assumed: OH 6498.729 air
converts to 6500.525, and `OH_GMOS_lines.dat` carries 6500.521 (0.004 A). Applying
an air/vacuum conversion here would fake a ~2.7 A offset at 10000 A — larger than
the error being looked for.

Result, on the three feige110 frames (123 OH lines in band, spacing 5.53 A):

```
frame 1   -0.60 A      frame 2   -0.60 A      frame 3   +0.20 A
tolerance  1.38 A  (25% of a line spacing)          -> AGREE
```

**The wavelength scale is confirmed against the sky to better than 1 A**, i.e.
well under the ~4 A a one-line misidentification would produce. Since the ThAr
catalogue, the template and the arc are all one closed loop, this is the only
number in the run that comes from outside it.

Read the *correlation* carefully: the peak value is only 0.14, which looks like a
failed match and is not. OH catalogue intensities do not predict observed line
strengths well, and the sky is faint in a 50-240 s standard exposure. What matters
is that the peak stands above the lag curve's own noise — 8.0 sigma on one frame,
14.4 sigma on the three stacked. `check_skylines_xcorr.py` reports that sigma
alongside the correlation for exactly this reason, and flags a peak under 4 sigma
as no measurement at all.

### The cube edge ring

Both cubes come out (2193, 15, 23) — a 31" x 20" Large-slicer field, the normal
geometry, and the same 15 x 23 grid as the known-good RL 2024-12-04 B cube. **RH4
did not inherit RH2's oversized-cube pathology** (553-1307 spaxels on one axis);
those RH2 configs were dithered, and RH4's frames share one pointing exactly.

feige110's cube does carry a spatial artifact: 18.1% of its signed flux sits in
the outermost row and column (rows 0/14 hold 10.7%/6.9%, columns 0/22 hold
11.5%/9.8%). The star itself is textbook — peak at spaxel (7,11), the grid centre
and the same spaxel as the RL cube, with 63% of the signed total inside **one**
spaxel. The ring is 4-11 spaxels from the peak, entirely outside the 3.4"
(2.5 spaxel) aperture, and §6 shows by direct comparison with feige34 that it does
not reach the throughput.

**This exposed a flaw in the guard inherited from `post_rh2.py`**, worth recording
because it cried wolf: that check measured flux within a fixed +/-8 spaxels of the
peak. On RH2's grids of hundreds of spaxels that is a core; on a 15 x 23 grid it
spans nearly the whole field, so it measured the edge ring and reported feige110
as a "smeared or doubled stack" at 77.6%. `run_rh4_throughput.py` now measures
concentration inside the **actual extraction aperture** and reports the edge ring
separately:

```
            old metric      new metric (2.5 spaxel aperture)   outer ring
feige110      77.6%   ->              76.9%                      18.1%   NOTE
feige34       99.2%   ->              97.9%                       0.1%   clean
```

The two metrics barely differ on feige110 — the point is not the number but that
the ring is now reported as its own quantity, so a real smeared stack (which would
drop the aperture number) stays distinguishable from an edge artifact (which does
not).

**A geometric threshold tuned on one grating's cube geometry is not a property of
the code** — the same lesson as the seed-span window in §2, in a different place.

### Warnings that are not problems

**96 `Right slit edge shift ... exceeds the maximum allowed of 0.0%`.** KCWI sets
`tweak_slits_maxfrac = 0.0` deliberately ("Make sure the full slit is used"), so
*any* computed shift trips it and the edge is left untweaked — which is the
intent. RH1's logs carry 704 of these and RH2's 238. Not RH4-specific, not a
defect.

**24 slits, 0 masked, widths 69.2-71.5 px.** No sign of the dead-column pathology
that invented a 25th slice on RH1's 2023-11-08.

### The coadd is slow here, and it is the machine, not the data

`pypeit_coadd_datacube` took roughly two hours per cube against the ~9 minutes
THROUGHPUT_PROCEDURE.md quotes. That was chased on the assumption it meant an
oversized grid, and it does not: the grids are the normal 15 x 23. The cause is
memory — the process holds ~2.4 GB and the 16 GB box was already swapping
(3.7 of 5.1 GB swap in use), so the per-slit rate collapsed from ~2 minutes to
~15 and recovered when pressure eased. RH4's cube has 2193 wavelength planes
against RH2's ~2090, so it is marginally larger, not structurally different.

Recorded so the next person does not go looking for a data defect. If it needs to
be faster, close other applications before starting rather than changing the
subpixel sampling, which would alter how flux is distributed between spaxels.

### The stray DARK config

`pypeit_setup` puts the night's one DARK (dispname `None`, cenwave 0.0) in its own
setup A with seven biases and nothing else. The incomplete-config guard inherited
from `run_rh2.py` skips it. No action needed; noted because a setup letter that
holds nothing is otherwise alarming.

### Dome flats were typed correctly

Unlike RH2's 2024-11-06, this night's three dome flats carry `FLSPECTR = 'on'` and
PypeIt typed them `illumflat,trace` unaided. No hand edit, so — unlike the RH2
run — `pypeit_setup` here is safe to re-run.

---

## 6. Results

```
cube                    star      cenw  airm   fit range (A)   interior mean   peak
feige110_2024-12-28_B   feige110  9950  1.20   9479.9-10179.8      32.1%    34.2% @ 9778 A
feige34_2024-12-28_B    feige34   9950  1.15   9480.2-10235.8      31.8%    34.5% @ 9757 A
```

Fit quality, on the RL/RH2 metrics (proportional windows, per RH2 §4.4):

```
                    feige110   feige34    RL order-15 reference
interior scatter     0.0260     0.0259           0.0300
worst edge error     0.083      0.079            0.117
```

Both are *better* than the RL reference, on configs spanning ~700-750 A against
RL's 3400 — so polyorder 15 is not overfitting here either.

Plot: [rh4_throughput_curves.png](../results/rh4_throughput_curves.png).

### The curve is shaped by the detector, not the blaze

```
  9500 A  27.2%      9800 A  34.2%  <- peak      10000 A  25.8%
  9600 A  31.4%      9900 A  31.4%               10100 A  18.8%
  9700 A  33.7%
```

RH2 found its throughput peak tracking the grating angle, 38 A from cenwave on
average, as a blazed grating should. **RH4's peak sits 172 A blueward of its
9950 A cenwave and the curve then falls steeply** — 34.2% to 18.8% in 300 A.
That is the silicon QE collapse approaching 1 micron, not the blaze: it is far
too steep and too red-sided to be a blaze function, and it is monotonic past the
peak. RH4 is the grating where the detector, not the grating, sets the envelope.

Read the extreme ends with the usual caution — 10100 A is 80 A from the fit
boundary, inside the ~150 A zone THROUGHPUT_PROCEDURE.md marks as good to only
~12%.

### Cross-check against RL — partial agreement, and the overlap is thin

The only external check this single night affords. It is weaker than it first
looks, and the reason is worth stating rather than hiding in a plot.

RL's composite runs to 9860 A, but **it is supported by two or more configs only
out to 9755 A** — beyond that it is a single config, and the `throughput_std`
column there (0.07%, then 0.01%) is the standard deviation of one number, not a
measurement of agreement. RH4's own reliable interior starts around 9655 A
(9480 + a quarter of its range). So the genuinely comparable window is about
**9655-9755 A — roughly 100 A**, not the 380 A the ranges suggest.

Inside it:

```
          RH4     RL    diff   RL n_configs
 9600 A  31.4%  32.7%   -1.3        2
 9650 A  32.8%  32.5%   +0.3        2
 9700 A  33.7%  31.9%   +1.8        2
 9750 A  34.2%  30.9%   +3.2        2
```

They cross at ~9650 A and diverge to +3.2 points by 9750. Both curves are moving
fast and in opposite directions there (RH4 rising to its peak, RL falling into
its red end), so some divergence is expected — but a 3-point gap is larger than
either measurement's internal scatter, and it is **not** resolved by this data
set. RL's red end is also its own weakest region: THROUGHPUT_PROCEDURE.md records
that a red-edge throughput "decline" on RL was previously traced to
partial-coverage bias, which is a reason to suspect the RL side of the gap rather
than the RH4 side, but that is a hypothesis, not a measurement.

**Conclusion: consistent at the ~1-point level where the overlap is best
(9650-9700 A), unresolved at the 3-point level by 9750 A.** Not a clean
confirmation, and it should not be quoted as one.

### The two standards agree — the check that does pass cleanly

This is the strong result of the run. feige110 and feige34 were observed 11 hours
apart, sit at different places in the IFU, and were reduced through separate sky
regions, separate cubes, separate extractions and separate sensfuncs. Their
throughputs:

```
          feige110   feige34    diff
 9550 A     28.8%     29.2%     -0.4
 9650 A     32.8%     33.0%     -0.1
 9750 A     34.2%     34.5%     -0.3
 9850 A     33.5%     33.8%     -0.3
 9950 A     28.6%     29.7%     -1.1
10050 A     22.2%     24.2%     -2.0

common interior 9655-10005 A:   32.1%  vs  32.6%   ->  0.6 points
```

Agreement is 0.1-0.4 points from 9550 to 9850 A, widening only past 9950 where
both curves are plunging down the QE cliff and small wavelength differences turn
into large throughput differences. **A different star reproducing the same
throughput tests the flux calibration itself**, which is exactly the cross-check
RH2 got from g191b2b, and it is the one independent confirmation available on a
single night.

It also settles the cube-edge question of §5 by measurement rather than by
argument: feige110's cube carries **18.1%** of its signed flux in the outermost
row/column and feige34's carries **0.1%**, yet the two throughputs agree to 0.6
points. The edge ring does not propagate into the extraction.

---

## 7. Still open

**The 3-point gap against RL at 9750 A** (§6). Not resolved by this data set. The
cheapest test is a second RH4 night at a different central wavelength: if RH4's
own curves at two grating angles agree with each other where they overlap, the
gap points at RL's red end rather than at RH4.

**A second RH4 night would also fix the gate** (§3). Everything RH4 knows about
its own wavelength scale comes from one arc, one template seeded from that arc,
and a catalogue chosen for the band. The sky check (§5) is the only outside
measurement, and it constrains a common-mode shift to ~1 A but says nothing about
a second night's reproducibility.

**Whether `ThArRH4_lines.dat`'s amplitude cut is near-optimal is untested.** 100
was chosen to put catalogue density near observed line density and it worked
first time (24/24 slices, RMS 0.290 px). Nobody has checked whether 50 or 200
does better. The `--amp-min` flag exists for that.

---

## Tools added for RH4

| file | purpose |
|---|---|
| `pypeit_test/build_rh4_linelist.py` | write `ThArRH4_lines.dat`; the piece RH1 and RH2 never needed (§1) |
| `RH4 pypeit run/seed_rh4.py` | retype ThAr as the arc, holy-grail seed; tolerates the expected `rc=1` |
| `pypeit_test/build_rh4_template.py` | build `keck_kcrm_RH4.fits`; single seed, no stitch (§3) |
| `pypeit_test/gate_rh4_template.py` | gate template **and** line list together; prints the self-consistency caveat |
| `RH4 pypeit run/run_rh4.py` | night driver; sky regions per star, one reduction pass per star (§4) |
| `RH4 pypeit run/run_rh4_throughput.py` | cubes split by star, extraction, trim, sensfunc |
| `pypeit_test/check_skylines_xcorr.py` | sky-line audit by cross-correlation, at any wavelength (§5) |

Changed and shared with the earlier runs:

| file | change |
|---|---|
| `pypeit_test/pick_seed_slits.py` | `rank()` and `health()` take `span=(MIN, MAX)`; CLI takes `--span`. Defaults unchanged, so RH1/RH2 behave as before (§2). Backup at `.bak`. |
| `pypeit_test/gate_rh2_template.py` | fixed a pre-existing crash: it unpacked 4 values from `health()`, which has returned 5 since the shifted-slice check was added. It now also reports shifted slices. |
| `pypeit_test/plot_rh1_throughput.py` | none needed — `--grating RH4` works as written. |
