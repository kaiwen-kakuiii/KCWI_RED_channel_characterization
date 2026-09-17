# RH2: what differed from the RL and RH1 runs

Only the deltas. The pipeline itself is [PYPEIT.md](PYPEIT.md) and
[THROUGHPUT_PROCEDURE.md](THROUGHPUT_PROCEDURE.md); the missing-template argument
is [RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) and is not repeated here. This file is
the problems RH2 raised that those do not cover, and what was done about each.

Written 2026-08-24, PypeIt 2.0.1, spectrograph `keck_kcrm`.

Scope after the drops below: **6 nights, 8 configs, 29 science frames** —
feige110 on five nights, g191b2b on 2023-09-16. Every config solved 24 of 24
slices.

---

## 1. Two templates, not one — the slicer splits the run

RH2 has the same hole RH1 does: `keck_kcwi.py::config_specific_par` assigns
`reid_arxiv` for RL, RM1, RM2 and RH3 only, so RH2 falls back to holy-grail. That
part is [RH1_WAVELENGTH.md](RH1_WAVELENGTH.md).

**What RH1 did not face: RH2's nights use two different slicers**, and the slicer
sets spectral resolution. A template built at one slicer holds as single features
the blends the other resolves, so the line inventories differ and identification
fails — exactly how RH1's Small-slicer setup died against a Large-slicer
template. So `build_rh2_template.py` **requires** `--decker` rather than
defaulting it; building one template across both would silently mix resolutions.

| template | slicer | shape | seeds |
|---|---|---|---|
| `keck_kcrm_RH2.fits` | Large | 2918 px, 6492.1-7554.2 A, 0.3657 A/px | 2, stitched at 7034.6 A |
| `keck_kcrm_RH2_medium.fits` | Medium | 2064 px, 6358.2-7159.0 A, 0.3895 A/px | **1, no stitch** |

### Does it need a stitch? Ask per slicer, not per grating

One RH2 exposure covers **~764 A per slice** (measured, 2023-11-17 B best slice:
6492.1-7255.8 A). Whether that is enough depends entirely on how the slicer's
nights are spread:

```
LARGE — three central wavelengths, union outruns one exposure
     6400      6600      6800      7000      7200      7400      7600 A
       |         |         |         |         |         |         |
6900 [=================================]
7100         [=================================]
7180             [=================================]
       |<--------------- union 1086 A ------------->|
       |<-- one exposure 764 A -->|                        -> STITCH NEEDED

MEDIUM — every night at the identical cenwave 6749.978
6750 [=================================]                   -> ONE SEED, NO STITCH
```

The Medium case is the good one and worth recognising when it appears: a
single-seed template is one exposure wide, which is the shape PypeIt ships, and
it has no stitch joint and no alignment ambiguity. Built with
`--max-seeds 1 --target 6343 7159`.

Coverage is the thing that must not fall short, so the Medium seed was chosen as
the **widest** slice (spat 1901, rms 0.296) rather than the lowest-RMS one (spat
1672, rms 0.083). That is deliberate: at 0.3895 A/px an rms of 0.296 px is
0.115 A against a line spacing of ~8.6 A, i.e. irrelevant to identification,
while 24 A of extra coverage is not.

### Holy-grail did better here than on RH1

RH1 got 1 slice of 24 and could not run the flat field at all. RH2:

| seed setup | slicer | cenwave | usable slices |
|---|---|---|---|
| 2023-11-17 B | Large | 6900 | **7 / 24** |
| 2023-07-16 A | Large | 7180 | 2 / 24 |
| 2024-11-07 B | Medium | 6750 | 4 / 24 |

Enough that seeding took one pass with no bootstrap round. Each run still ends
`rc=1` — the unsolved slices kill the flat field — which is expected, because
WaveCalib is written before the flat field. `seed_rh2.py` treats non-zero as
normal and only checks that a WaveCalib appeared.

### The gates

Both passed 24/24 with `run_pypeit -c` returning 0:

| template | gated on | slices | median rms |
|---|---|---|---|
| Large | 2023-11-17 B | 24/24 | 0.058 px |
| Medium | **2024-10-29 B** | 24/24 | 0.269 px |

**The Medium gate is a stronger test than RH1's was.** RH1 gated on a setup its
seeds came from; the Medium template was gated on a *different night* from its
seed, so it demonstrates the template generalises rather than reproducing its own
arc. Worth copying for any future grating.

Medium's rms is 4.6x Large's, entirely from the deliberately wide seed above.
0.269 px = 0.105 A, and the RL work already recorded that a 0.20 A wavelength
error is irrelevant to a throughput continuum. Confirmed downstream: the Medium
nights' sky models and throughput scatter are as good as the Large nights'. It
was **not** rebuilt, because swapping a passing template without re-gating is how
hours get lost.

### Seed agreement, for the record

`check_seed_agreement.py` on the two Large seeds, over their 442 A overlap:

```
cross-correlation peak 0.923 at +0.304 A   (0.04 line spacings)
zero-shift correlation 0.866
AGREE: 0.30 A against a 2.15 A tolerance (25% of a line spacing)
```

Much tighter than RH1's (+1.13 A, 0.27 spacings, zero-shift correlation 0.190).
The RH2 seeds were nearly aligned before any shift.

---

## 2. Four data problems the RL and RH1 nights never showed

Each was caught by inspection *before* spending a reduction on it, and each is
now guarded in code rather than in memory.

### 2.1 Dome flats whose lamp keyword lies (2024-11-06)

PypeIt typed that night's three dome flats `None`, leaving the night with **zero
trace frames** — and it is the night's only config, so the night was unreducible.

`check_frame_type` requires `dome_noarc` lamps on, and the header says
`FLSPECTR = 'off'`. The pixels say otherwise:

```
2024-11-06 domeflats   24317, 24263, 24267 ADU      FLSPECTR = off
2024-10-29 domeflats   24575, 24038, 24068 ADU      FLSPECTR = on    (known good)
```

Same exposure time, same illumination level. A bookkeeping error, not dark
flats — so the three rows were retyped `illumflat,trace` by hand. Vindicated: the
night then solved 24/24 slices and its extracted counts (BOX_COUNTS median 605)
land beside the other Medium nights (612, 619).

**Rule:** when a frametype looks wrong, measure the counts before believing the
keyword.

### 2.2 A zero-second frame, and a 5-second one

2024-11-07 carried a `W1` frame with `ELAPTIME = 0.0` (aborted, and not a
standard). 2023-11-07 carried a 5 s feige110 frame that sorted **first** among
its science rows — so it would have been the frame handed to the object finder.

Both are now handled in `run_rh2.py` rather than by hand:

```python
MIN_SCI_EXPTIME = 10.0                       # below this is an acquisition/test frame
best = max(sci, key=lambda f: exptime.get(f, 0.0))   # finder gets the LONGEST exposure
```

`sci[0]` was the RH1 driver's choice and carried no information about frame
quality. The finder now runs on 60-600 s frames throughout.

### 2.3 Incomplete configs from stray frames

`pypeit_setup` groups on (dispname, decker, binning, cenwave), so a frame whose
header disagrees with its night lands in a setup of its own with no calibrations.
2023-11-17 setup D is one science frame, `decker: unknown`, seven biases, nothing
else. Running it wastes an hour to reach a confusing error, so the driver skips
any setup missing arc/trace/pixelflat and says so.

### 2.4 A night with no star in it (2023-07-16)

`find_object_regions.py` returned `user_regions = :0,100:`. That string is
**degenerate**, and the failure is silent: PypeIt's `read_userregions` parses
`:0` as `[0,0]` and `100:` as `[res-1,res]`, so **no sky pixels are defined at
all** — and the reduction still completes. `check_skysub.py` then reports `nan`
sky-zone residuals with a meaningless "sky ok" verdict on all 24 slices.

The frame was checked by eye: no star, target unknown. Night dropped.

The finder now refuses to emit a degenerate region:

```
NO DETECTION: object would span 0%-100% of the slice (100% wide, degenerate at >=85%).
  No sky region can be defined, so no user_regions is suggested.
  ...
```

and `run_rh2.py` stops that setup, logs it to `NEEDS_DECISION.txt`, and continues
with the rest rather than reducing it. Threshold margin: every real night runs
26-58% object width (RH1's widest was 58%), against the 85% gate.

**A `:0,100:` from the finder means it found nothing, not that it misfired.**

### Also dropped

**2024-11-03** — poor weather, Kaiwen's call. Raw files remain in `fits/lev0/`;
`by_night` is only a hardlink view, so a deleted night looks exactly like an
accident. Do not restore either night.

---

## 3. Driver changes against `run_rh1.py`

`RH2 pypeit run/run_rh2.py` differs in four ways, all forced by the above:

1. **Template chosen per setup from its own `decker`**, not one global constant.
2. **Incomplete configs skipped by inspection** (§2.3).
3. **`MIN_SCI_EXPTIME` and longest-exposure frame selection** (§2.2).
4. **`pypeit_setup` runs only when no `.pypeit` exists.** Two files carry hand
   edits (§2.1, §2.2) that regenerating would silently revert — taking a night's
   only complete config with it.

Retained from RH1 unchanged: the `"reid_arxiv" not in txt` test before patching
(so a setup already pointed at the *other* slicer's template is not repatched and
its calibrations deleted), and the rule that calibrations are keyed by frame, not
by parameters.

---

## 4. Post-processing: what THROUGHPUT_PROCEDURE.md does not cover

### 4.1 `align = True` is a red herring — do not repeat this

2023-11-07 B is the one dithered config (a 2" nod, DECOFF 0 and -2). The doc says
to set `align = True` when frames are offset. Its cubes came out with one spatial
axis tens of times too large, and alignment looked like the obvious cause.

It is not. Rebuilding with `align` removed gave a **byte-identical cube geometry**
and BOX_COUNTS within 0.5% (637 vs 634). The header WCS already places a small
nod correctly.

**`align = True` is for when the WCS cannot be trusted, not merely when frames are
offset.** `post_rh2.py` writes no `align` line for any config.

### 4.2 The oversized cubes, still unexplained

Three feige110 Large configs produced cubes with one spatial axis 553-1307
spaxels where ~15 is correct:

| config | star | slicer | grid | field |
|---|---|---|---|---|
| 2023-11-07 B | feige110 | Large | 24 x 1289 | 33" x 1750" |
| 2023-11-17 B | feige110 | Large | 553 x 23 | 751" x 31" |
| 2023-11-17 C | feige110 | Large | 1307 x 23 | 1775" x 31" |
| **2023-09-16 B/C** | **g191b2b** | **Large** | **23 x 15** | **31" x 20"** |
| all Medium configs | feige110 | Medium | 30 x 23 | 20" x 16" |
| RL 2024-12-04 B | feige110 | Large | 15 x 23 | 20" x 31" |

(Grid and field are both `NAXIS2 x NAXIS1`. The g191b2b and RL Large cubes are
the same 20"x31" field with the axes swapped — the rotator angle decides which
sky axis the slices run along, and it carries no other meaning here.)

So it is **not** "RH2 + Large" — the g191b2b Large configs are normal. Cause
unknown; it is specific to those three nights.

**It is benign for throughput, and that was measured, not assumed:** on
2023-11-17 B the star peaks at spaxel (7,11) — the same spaxel as the good RL
Large cube — and rows 0-15 hold **99.82%** of the flux. The other 537 rows carry
0.18%.

`post_rh2.py` therefore **warns** on poor flux concentration (<90% within +/-8
spaxels of the peak) but **rejects only** a degree-scale `CDELT`, which is the
documented single-frame failure that inflates counts by ~1e21. An earlier version
rejected on field size and would have thrown away every valid RH2 cube.

Caveat on that concentration metric: a sky-subtracted cube contains negative
pixels and the denominator is the signed total, so values are not monotonic in
aperture and can exceed 100% (one config reads 185%). Read it as "compact,
comparable to RL", not as aperture photometry. Where it looked alarming
(2023-09-16 C, 185%) `check_skysub.py` confirmed the sky model is sound —
STD_CHIS 1.019-1.128, the lowest MED_CHIS of any RH2 config.

### 4.3 THROUGHPUT is not a column

In the sens file, `SENS` is a table but throughput is **`ImageHDU` 5**, beside
`WAVE` (3) and `ZEROPOINT` (4). Reading `d["THROUGHPUT"]` off the `SENS` table
finds nothing and fails silently.

### 4.4 polyorder 15 on a ~700 A range

The RL tuning chose order 15 across ~3400 A. An RH2 config's fit range is
~700 A — one degree of freedom per ~46 A against RL's ~227 A. It holds anyway:

```
interior scatter (mag)   RH2  0.016 - 0.051      RL order-15 reference  0.0300
worst edge error (mag)   RH2  0.090 - 0.204      RL order-15 reference  0.117
fit pixels per config    3493 - 3851  against 16 coefficients
```

Nowhere near interpolation, so the low values are real. **Do not reuse RL's metric
windows on RH2**: "interior = >300 A from each end" leaves 2800 A on RL and
~90 A here, and gave one config a meaningless 0.0042 mag. Use proportional
windows (middle 50%) instead.

### 4.5 Do not compare RH2 configs at a fixed wavelength

RH2's throughput peak **tracks the grating angle**; RL's does not:

```
mean |peak wavelength - cenwave|     RH2    38 A        RL   1343 A
```

Expected for a blazed high-dispersion grating — rotating it moves the blaze peak
with it — while RL spreads 3400 A and washes the structure out. Each RH2 curve
legitimately bows to its own centre, so a fixed-wavelength comparison across
configs measures the blaze, not the instrument. `plot_rh1_throughput.py` already
colours by central wavelength for this reason; it now takes `--grating RH2`.

---

## 5. Results

```
cube                    cenw  airm   range (A)      interior mean
feige110_2024-11-06_B   6750  1.11  6404.6-7127.8      24.4%
feige110_2024-11-07_B   6750  1.29  6407.2-7132.6      24.9%
feige110_2024-10-29_B   6750  1.38  6407.5-7132.0      25.5%
feige110_2023-11-17_B   6900  1.17  6543.4-7208.9      30.0%
g191b2b_2023-09-16_C    6950  1.19  6596.7-7296.5      30.6%
feige110_2023-11-17_C   7100  1.15  6759.1-7438.3      32.1%
g191b2b_2023-09-16_B    7100  1.20  6750.6-7439.7      31.0%
feige110_2023-11-07_B   7100  1.24  6749.2-7438.5      29.8%
```

Plot: `rh2_throughput_curves.png`.

**Cross-check that passes:** three configs at cenwave 7100, two different
standards, three nights over 15 months — 29.8%, 31.0%, 32.1%, with g191b2b
between the two feige110 values. A different star reproducing the same throughput
tests the flux calibration, and it holds.

The g191b2b curves are the least precise of the eight (interior scatter 0.0513
and 0.0358 mag, edge error 0.204 and 0.142) — 90 s exposures against feige110's
240-600 s.

### Open: the 6750 offset is confounded

Throughput rises with central wavelength, ~25% at 6750 to ~31% at 7100. It is
**not established** whether that is the grating or the slicer, because in this
data set the two cannot be separated: every 6750 config is Medium and every
6900/6950/7100 config is Large. Comparing at each config's own blaze peak leaves
the gap essentially unchanged (Large 31.7% +/- 1.1, Medium 26.2% +/- 0.3), so the
blaze does not explain it either.

**The test exists and was not run:** KOA holds g191b2b RH2 frames on
**2024-12-25** at cenwave 6700 and 7000, **Large slicer** — a Large config 50 A
from the Medium configs' 6750. If its throughput lands near 26%, central
wavelength drives the trend; near 31%, something is slicer-dependent. That night
was deliberately not fetched.

---

## Tools added for RH2

| file | purpose |
|---|---|
| `RH2 pypeit run/seed_rh2.py` | holy-grail seed runs, per slicer; tolerates the expected `rc=1` |
| `pypeit_test/build_rh2_template.py` | stitch a template; `--decker` is **required** |
| `pypeit_test/gate_rh2_template.py` | run one setup against a template and count solved slices |
| `RH2 pypeit run/run_rh2.py` | night-by-night driver (§3) |
| `RH2 pypeit run/post_rh2.py` | writes each `.coadd3d`, coadds, extracts; cube guards (§4.2) |
| `RH2 pypeit run/sens_rh2.py` | trim + sensfunc, strictly serial (telluric cache race) |
| `pypeit_test/plot_sens_overlay.py` | all-config zeropoint/throughput overlay; also works on the RL run |

Changed and shared with the RH1 run:

| file | change |
|---|---|
| `pypeit_test/find_object_regions.py` | `NO DETECTION` guard at >=85% object width (§2.4). Fires only in the degenerate case; RH1's widest real night was 58%, so RH1 is unaffected. Backup at `.bak`. |
| `pypeit_test/plot_rh1_throughput.py` | takes `--grating RH2`; RH1 remains the default. Also fixed the RL reference band's label being pinned at 8600 A, off the end of RH2's range. |
