# What RH1 needs that RL did not

Companion to [PYPEIT.md](PYPEIT.md) and [THROUGHPUT_PROCEDURE.md](THROUGHPUT_PROCEDURE.md),
which describe the standard path. **Everything in those two still applies.** This
file is only the delta: the six things that broke on RH1 and did not break on
RL, the one extra pipeline step, and the checks added because the existing ones
passed defective data.

Written 2026-08-25, against PypeIt 2.0.1.

Scope: 9 nights, 10 configurations, 44 standard-star frames. The RH1 standards
are **g191b2b** (36 frames, 7 nights) and **feige34** (8 frames, 2 nights) —
not the feige110 of the RL run. Both are calspec standards within 2" of their
catalogue positions, so `pypeit_sensfunc` matches them on coordinates unaided.

---



## Vocabulary

Terms this file leans on, defined once. Terms shared with the RL documents are
not repeated.


| term                        | meaning                                                                                                                                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **central wavelength**      | wavelength at the middle of the detector, set by the grating angle (header `RCWAVE`). Changing it rotates the grating                                                                                         |
| **blaze**                   | the grating's efficiency peak. It moves with the grating angle, so each central wavelength has its own efficiency curve                                                                                       |
| **template** (`reid_arxiv`) | a stored arc spectrum *with* its wavelength solution, used to decide which catalogue line each observed feature is                                                                                            |
| **holy-grail**              | PypeIt's fallback that identifies lines from scratch, with no template                                                                                                                                        |
| **slicer / decker**         | the image slicer: Small, Medium or Large. Its width sets the spectral resolution                                                                                                                              |
| **spaxel**                  | one spatial pixel of the reconstructed datacube. **Its angular size depends on slicer and binning** — 1.358" here at Large/2x2, 0.339" at Small/1x1 — so a radius in spaxels is not comparable between nights |
| **enclosed fraction**       | share of the star's whitelight flux falling inside the extraction aperture. Quoted at a stated radius in **arcsec**, never spaxels                                                                            |
| **visit**                   | a block of frames taken close together. Five RH1 setups were observed in two visits hours apart                                                                                                               |
| **shifted slice**           | a slice whose wavelength solution is displaced bodily from its neighbours' — normal width, normal dispersion, wrong wavelengths                                                                               |


---



## The six problems, at a glance


| #   | symptom                                   | root cause                             | fix                    |
| --- | ----------------------------------------- | -------------------------------------- | ---------------------- |
| 1   | calibration crashes in the flat field     | PypeIt ships no RH1 template           | build one              |
| 2   | Small-slicer night solves 0 of 25 slices  | template built at the wrong resolution | a second template      |
| 3   | 4 slices land 201-337 A off               | template wider than one exposure       | one-exposure templates |
| 4   | spurious 25th slice, then NaN in the cube | a dead detector column                 | `exclude_regions`      |
| 5   | Small-slicer night reads 2.4 pts low      | extraction aperture defined in spaxels | `--boxcar 3.4`         |
| 6   | throughput scatter 3x the RL result       | fewer photons genuinely arrived        | not a reduction fault  |


Problems 1-3 are all the same underlying subject — the template — approached
from three directions. Problem 4 is unrelated. Problem 6 is a result, not a bug.

**Problems 2 and 5 are the same mistake in two places**: a quantity expressed in
detector or grid units, compared across configurations where the unit is not the
same size. Arc line width in pixels instead of Angstroms; aperture radius in
spaxels instead of arcsec. Both cost a day. When a number carries a unit that
belongs to the instrument rather than to the sky, convert before comparing.

---



## 1. PypeIt ships no template for RH1

`keck_kcwi.py::config_specific_par` assigns a template per grating, and the list
has a hole:


| grating           | template                            |
| ----------------- | ----------------------------------- |
| BL, BM, BH2, BH3  | `keck_kcwi_*.fits`                  |
| RL, RM1, RM2, RH3 | `keck_kcrm_*.fits`                  |
| **RH1, RH2, RH4** | **none — falls back to holy-grail** |


Not a local problem; no `keck_kcrm_RH1.fits` exists on PypeIt's `develop` branch
either. **RL worked precisely because RL is on the supported list**, which is why
nothing in the RL documents prepares you for this.

Holy-grail solved **1 slice of 24**, and the 23 unsolved slices then crashed the
flat field:

```
ValueError: zero-size array to reduction operation minimum which has no identity
  flatfield.py:2220  mnmx_wv[slit_idx, 0] = np.min(waveimg[onslit_init])
```

Its failures are not honest. Two slices *reported* a better RMS than the one
worth keeping — 0.129 and 0.192 — through six points spanning 10^5 A. A returned
fit object is not a solution.

**Fix.** Bootstrap: run holy-grail once, keep the single slice that survives
vetting, build a template from it, and that template then solves all 24. See
[RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) for the full account.

Gate on the night that crashed:


|            | slices solved | median RMS   | outcome |
| ---------- | ------------- | ------------ | ------- |
| holy-grail | 1 / 24        | 0.346 px     | crashed |
| template   | **24 / 24**   | **0.123 px** | clean   |


The clearest demonstration is 2023-10-15: holy-grail gave it 1 slice of 24 at
RMS 0.499, and the *identical frames* against a template gave 24 of 24 with a
best slice at 0.064. Same data; only the line identification changed.

---



## 2. One template is not enough — the slicer sets the resolution

2023-11-08 is the only **Small**-slicer setup, and the Large-slicer template
solved **0 of its 25 slices**.

**Compare arc line widths in Angstroms, never in pixels.** This is the single
most useful lesson here, and getting it wrong cost most of a day. In pixels the
failing setup and a known-good RL one look identical — both ~2.7 px. The pixel
is a different physical size in each case:


| spectrum                               | FWHM px | A/px   | **FWHM (A)** | template/data |
| -------------------------------------- | ------- | ------ | ------------ | ------------- |
| shipped `keck_kcrm_RL`                 | 4.00    | 1.8356 | 7.34         | —             |
| RL 2023-12-17, Small, 1x1 — **works**  | 6.50    | 0.9181 | 5.97         | **1.23x**     |
| our RH1 template (Large seeds)         | 5.00    | 0.2973 | 1.42         | —             |
| RH1 2023-11-17, Large, 2x2 — works     | 5.00    | 0.2841 | 1.42         | 1.00x         |
| RH1 2023-11-08, Small, 1x1 — **fails** | 4.50    | 0.1531 | **0.69**     | **2.06x**     |


The RL Small-slicer night works because RL is low-dispersion: 5.97 A lines
against a 7.34 A template, near-matched. RH1 at the Small slicer reaches 0.69 A
against a 1.42 A template, so the data **resolves blends the template holds as
single features** — the two line inventories are not the same list.

**Fix.** A second template, `keck_kcrm_RH1_small.fits`, seeded by holy-grail from
that setup's own arcs. Holy-grail does markedly better at the Small slicer
(13 of 25 slices, against 1 of 24 at Large), because cleaner line separation is
exactly what its pattern matching needs.

**So: one template per grating AND per slicer.** A KCRM programme using two
slicers needs two templates.

---



## 3. A template must be one exposure wide

This one is counter-intuitive and cost two wrong diagnoses.

One RH1 exposure covers ~630 A while the nine nights span 5800-7000 A, so the
first template was stitched to 1103 A — 1.8 exposures — to cover everything.
That width is a defect, because of how `full_template` works:

```python
i0    = npad // 2 + int(shift_cc)      # one cross-correlation picks the position
mspec = temp_spec[i0 : i0 + nspec]     # a window as wide as YOUR data
mwv   = temp_wv  [i0 : i0 + nspec]     # ...and the wavelengths that go with it
```

The template is not used whole. A window the width of your data is **cut out of
it**, at a position chosen by a single cross-correlation.

```
template   5796 |==============================| 6899    3720 px = 1.8 exposures
data                6216 |==============| 6820          2064 px

freedom in where to cut:  1103 - 604 = 499 A  =  ~1725 px
```

Two consequences. Half the template has no counterpart in the data, which
dilutes the correlation to cc~0.24 for *every* slice; and 1725 px of freedom
with that flat a correlation surface means noise picks the winner.

Measured on 2023-10-17 C, where 22 slices chose one alignment and 2 chose another
~1090 px away (~315 A):


| slice                | shift chosen        | cc              |
| -------------------- | ------------------- | --------------- |
| 516, 884, 1163, 1236 | +455 to +547 px     | 0.239-0.248     |
| **1310, 1458**       | **-541 to -549 px** | **0.272-0.274** |


**The wrong choices scored higher.** Nothing was malfunctioning; the correlation
genuinely preferred the wrong window.

Four slices across two setups, both at central wavelength 6520 — the setting
furthest from the template's blue start:

```
2023-10-17 C   spat 1310 (-331 A),  spat 1458 (-337 A)
2023-10-20 B   spat  812 (-337 A),  spat 1312 (-201 A)
```

**Why it matters downstream.** `measure_trim --spec2d` keeps only the interval
covered by *all* slices, so two outliers cut the usable range from ~555 A to
262 A. And ~2.3% of the star's light lands 335 A away, biasing those cubes low.

**Fix.** Build a template one exposure wide, aimed at the setup being calibrated:

```bash
python pypeit_test/build_rh1_template.py \
    --only "<night>/pypeit_run/<setup>" --max-seeds 1 --target 6216 6820 \
    --out pypeit_test/keck_kcrm_RH1_6520.fits
```

Before and after, on the same frames:


|                            | wide template | one-exposure template |
| -------------------------- | ------------- | --------------------- |
| shift spread across slices | ~1090 px      | **118 px**            |
| slice agreement            | 2 shifted     | **24 / 24**           |
| median RMS                 | 0.124 px      | **0.102 px**          |
| trim                       | 1243, 1124    | **126, 126**          |


**This is why every template PypeIt ships is exactly one exposure wide** — 2064 px
for a 2x2 KCRM read. Stitching solves coverage and buys an alignment problem in
exchange. Use narrow templates per central wavelength, not one wide one.

---



## 4. A dead detector column invents a slice

Unrelated to wavelengths, and specific to 2023-11-08.

Edge tracing found **25 slits where KCRM has 24**. The extra one was 4 px wide,
sitting between two real slices, and the CONTBARS frame found **zero** alignment
bars in it, which aborts calibration.

The cause, read off the dome flat:

```
counts by spectral band, spatial 650-655

rows          650    651    652    653    654    655
1800-2400    8911   8911   8892   8872   8842   8762
3000-3600    7236   7243   7222   7245   7192   7112
3600-4128    5923   5893   5898      6   5845   5633   <- column 653 dies
```

**Column 653 goes dead above spectral row ~3480** — 6 counts against ~5900 on
either side. PypeIt finds slit edges from the spatial gradient, and a drop from
5900 to 6 is a far sharper gradient than any real slit edge, so the tracer
registers an edge there and pairs it with the real edge at ~658.

**PypeIt's bad-pixel mask does not flag this column** — `spectrograph.bpm()`
returns zero flagged pixels in column 653, all 4128 rows.

Rejecting the sliver with `length_range` is not enough: it leaves the real slice
truncated at 130.4 px instead of 139.8, and a truncated slice makes the
astrometric transform degenerate, so `AlignmentSplines` divides by zero, the
slit-length spline is NaN for all 4128 rows, and the cube dies at
`create_wcs: int(NaN)`.

**Fix.** Exclude the column from slit tracing:

```ini
[calibrations]
    [[slitedges]]
        exclude_regions = 1:652:655,
```

The **trailing comma is required** — the parameter is a list, and without it
PypeIt iterates the string character by character and fails with
`not enough values to unpack (expected 3, got 1)`.

Result: **24 slices, every one 139.4-140.0 px**, no truncation, cube builds.

The range is deliberately narrow: it covers the dead column with a margin but
stops short of the real edge at ~658, which still has to be found.

**The star was never on this column** — 1678 px away, 0.000% of its flux — so
this costs no photometry. Its entire effect was structural.

---



## 5. The extraction aperture is measured in spaxels, not arcsec

Same shape of error as problem 2, in a different quantity, and it survived
undetected until every other correction had been made.

`pypeit_extract_datacube` sizes its boxcar from a 2D Gaussian fit to the
whitelight image, `core/datacube.py`:

```python
if boxcar_radius is None:
    nsig = 4
    wid = nsig * max(popt[3], popt[4])     # sigma from the fit -- in SPAXELS
else:
    wid = boxcar_radius / np.sqrt(arcsecSQ)
```

**A spaxel is not the same angle on every night.** PypeIt sizes the cube grid
from the slicer and binning:


| configuration | spaxel     | nights          |
| ------------- | ---------- | --------------- |
| Large, 2x2    | **1.358"** | 8 of the 9      |
| Small, 1x1    | **0.339"** | 2023-11-08 only |


Four sigma in spaxels ought to be scale-invariant — sigma in spaxels times the
scale is sigma in arcsec, and the grid cancels. Three things break the
cancellation:

1. **The fit floors sigma at 0.5 spaxel.** `bounds = ([0, 0, 0, 0.5, 0.5, ...])`.
  On the coarse grid that floor is 0.68", giving a 2.72" aperture — four of the
   fifteen cubes sat exactly on it. On the fine grid the same floor is 4x smaller
   and never binds.
2. **The coarse grid undersamples.** The star spans about one spaxel in a 23x27
  image, so the fitted sigma absorbs the spaxel itself and comes out broad.
3. **The profile is not Gaussian.** Four sigma of a true Gaussian encloses
  99.97%; measured enclosure ranged 0.80-0.99. Real flux sits in wings the core
   fit never sees.

Together these mean the Large-slicer nights were **right by accident** — their
inflated sigma bought a 2.7-3.9" aperture that happened to swallow the wings,
while the Small-slicer night fitted sigma correctly and was punished for it:


|                  | aperture used | flux enclosed | @6600 A |
| ---------------- | ------------- | ------------- | ------- |
| 2023-11-08 v1    | 1.81"         | 0.82          | 27.0%   |
| 2023-11-08 v2    | 2.17"         | 0.80          | 27.2%   |
| every other cube | 2.72-3.85"    | 0.94-0.99     | ~29.5%  |


**Fix.** `-b` is already in arcsec, so the aperture only has to be stated:

```bash
python "RH1 pypeit run/run_rh1_throughput.py" --boxcar 3.4 --skip-coadd
```

`--skip-coadd` reuses the existing datacubes; only the extraction changes, which
turns a multi-hour rebuild into minutes. After the fix all 15 cubes extract at
3.40" and 2023-11-08 reads **29.3% / 28.4%**, inside the cluster.

**Check the aperture fits the field before choosing it**, from the illuminated
pixels rather than the nominal slicer size — a radius that runs off the field
silently returns whatever is there. At 3.4" every RH1 cube is 100% on-field, the
tightest being 2023-11-08 with the star 7.8" from the nearest edge.

Two nights barely moved and are noted in *Still open*: **2023-10-20 B v2** and
**2023-11-08 A v2** enclose only 0.81 and 0.84 even at 3.4", where their own
siblings reach 0.95.

---



## 6. The one extra pipeline step: split cubes by visit

Not a bug — a step the RL procedure does not need, because the RL configurations
were all tight single sequences.

Five of the ten RH1 setups were observed in **two visits four to five hours
apart**, spanning airmass 1.21-1.63:

```
2023-10-15 B   5 frames, 290 min, airmass 1.24-1.60, TWO pointings
2023-10-16 B   6 frames, 297 min, airmass 1.25-1.61
2023-10-17 C   6 frames, 308 min, airmass 1.28-1.57
2023-10-20 B   6 frames, 245 min, airmass 1.21-1.51
2023-11-08 A   5 frames, 273 min, airmass 1.21-1.63
```

PypeIt applies the extinction correction **once per cube from a single header
airmass**, and disables the per-frame correction outright. Its own source says so:

> *"This could be wrong when combining multiple standard star exposures if the
> airmass of the standard star exposures is significantly different... the
> standard star exposures are assumed to have similar airmasses."*
> — `coadd3d.py`, where `extcorr_sort = 1.0` sits behind `if False:`

Measured against PypeIt's Mauna Kea extinction curve at 6200 A, the bluest RH1
coverage and so the worst case:


| airmass spread | flux error |
| -------------- | ---------- |
| 0.05           | 0.39%      |
| **0.10**       | **0.78%**  |
| 0.20           | 1.57%      |
| 0.36           | 2.84%      |


So `run_rh1_throughput.py` groups frames into cubes by **airmass**, not by time:

```bash
python "RH1 pypeit run/run_rh1_throughput.py" [--max-dam 0.10]
```

Airmass is the criterion because it is the quantity that sets the correction.
Time is only a proxy — frames minutes apart near the horizon can differ more in
airmass than frames an hour apart near transit. **Pointing is a hard split
regardless**: 2023-10-15's two blocks sit 10 detector pixels apart (`ra_off` 3.3
versus 0.0), which `combine = True` would smear rather than stack.

Result: 15 cubes from 10 setups, each spanning <=0.03 in airmass.

Note also `coadd3d.py` line 1225, `self.weights[ff] = 1.0` — **frames are
combined with equal weight regardless of exposure time.** Three RH1 setups open
with a short acquisition frame (20-30 s before 60-90 s). Checked and harmless
here: those frames sit within 1.6-3.1% of the long ones in counts per second,
with identical centroids. Unbiased, just noisier than optimal.

---



## 7. Checks added, and why each exists

Every one of these was written because an existing check passed defective data.


| check                                 | catches                                      | why it was needed                                                                  |
| ------------------------------------- | -------------------------------------------- | ---------------------------------------------------------------------------------- |
| `pick_seed_slits.py`                  | fits that look fine by RMS                   | holy-grail returned RMS 0.129 through 6 points over 10^5 A                         |
| `check_slice_wavelength_agreement.py` | a slice shifted from its neighbours          | a shifted slice keeps normal width, dispersion and RMS — only neighbours expose it |
| `check_seed_agreement.py`             | two solutions disagreeing where they overlap | a low RMS only scores a fit against the lines it chose                             |
| `check_skyline_wavelengths.py`        | a common-mode wavelength error               | the arc cannot audit itself; sky lines never entered the calibration               |
| `gate_rh1_template.py`                | a template unfit to run nights through       | "it converged" is not "it is right"                                                |




### Verifying a wavelength scale against the sky

The concern: the template was seeded by holy-grail and every slice solved against
it, so a **common-mode** error would leave every fit tight and every slice
agreeing. Night-sky emission lines can see what the arc cannot — they are in the
science exposure, their wavelengths are fixed by atomic physics, and nothing in
the reduction ever fit to them.

```
setup            lines   mean offset   scatter
2023-09-16         3       +0.63 A      0.04 A
2023-10-20         4       +0.50 A      0.01 A
2023-11-08         5       +0.41 A      0.02 A   ([OI] only)
```

A misidentification by one arc line would show as **~4.2 A**. These sit 6-8x
below that. The residual ~0.5 A is arc-to-science flexure, which PypeIt leaves
uncorrected (`[flexure] spec_method = skip`), and it is *rigid* — the 0.01-0.04 A
scatter proves the shape of the solution is right to a hundredth of an Angstrom.

**Use** `[OI]` **singlets, not OH.** OH sky features are molecular doublets: they
merge at the Large slicer's 1.42 A and a single catalogue value describes them,
but at the Small slicer's 0.69 A they resolve and the centroid lands on one
component. Measured on 2023-11-08: OH gave 1.01 A scatter where `[OI]` alone
gave **0.02 A**. The checker now bases its verdict on singlets and warns when
only blends are in range.

### Reading a health metric correctly

`pick_seed_slits.py` reports two numbers, and they answer different questions.
Seed criteria are deliberately strict, because a template seed should be
pristine; applying them as a health check calls a good setup a failure. On
2023-10-16, four slices came in at RMS 0.58-0.69 px against ~0.10 for the rest —
rejected as seeds, while PypeIt flagged none and reduced all 24. At 0.289 A/px
that is a 0.20 A error, irrelevant to a throughput continuum.

Similarly, `check_skysub.py` returns **SKY MODEL SUSPECT** on essentially every
standard-star frame, RH1 and RL alike — its 5%-of-sky threshold is routinely
crossed on short exposures where the sky is faint. **Read STD_CHIS, not the
verdict.** RH1 measures 0.95-1.03 against the RL control's 0.955-0.976; the RL
skysub study measured 3.235 for a genuinely corrupted model.

---



## 8. What the numbers say

All 15 cubes extracted through a fixed 3.4" aperture, two nights excluded as
non-photometric (below):

```
@6200 A   31.11% +- 0.56%   (n=2)
@6400 A   29.58% +- 1.46%   (n=8)
@6600 A   29.58% +- 1.52%   (n=11)
own fit centre  31.16% +- 1.84%   (n=13 cubes)
RL, for reference           ~29.5% +- 1%
```

The central value lands on the RL result. **The scatter is 1.5x worse**, and
that is the finding worth carrying forward.

### Two nights are non-photometric and are excluded

Four cubes cover ~6200 A and spanned 21.6-31.7% — a 10-point spread that no
blaze argument can absorb. Comparing **within one star**, so the expected flux
is identical, over 6050-6350 A with extinction removed:


| star    | night          | counts/s | throughput |
| ------- | -------------- | -------- | ---------- |
| g191b2b | 2023-09-16     | 932      | 30.2%      |
| g191b2b | **2023-11-07** | 793      | **25.3%**  |
| feige34 | 2023-11-17     | 1660     | 29.6%      |
| feige34 | **2024-03-13** | 1170     | **20.8%**  |


Throughput tracks counts, which for a fixed star at a fixed wavelength is close
to tautological. What it does establish is narrower and still useful: **the
deficit is already present in the extracted counts**, so sensfunc, the zeropoint
polynomial and the trim interval did not cause it.

Three further checks locate it outside the instrument:

- **Configuration is identical.** 2023-11-07 and 2023-09-16 match to seven
decimal places — `RCWAVE` 6199.9794922, `RGRANGLE` 104.355423, `RARTANG`
94.3950424, both Large/2x2/L2U2. Only `ROTPOSN` differs, which is sky
orientation. Nothing instrumental can explain it.
- **The loss is achromatic.** The feige34 pair differs by 10 A in central
wavelength, so it was tested rather than assumed: the ratio of the two curves
is flat at **0.698 +- 0.010** with a 4.9% tilt across 280 A. A moved blaze
would bend the ratio; a constant factor is light loss.
- **Aperture explains only part.** Corrected to total flux the two still sit 14%
and 22% low.

Excluding them changes @6200 A from 27.56% +- 3.94% to **31.11% +- 0.56%** —
from the noisiest point in the set to the tightest, on two different stars. The
blue end was never anomalous; the sample was. The exclusion lives in
`DROP` in `rh1_throughput_table.py`, which `plot_rh1_throughput.py` imports so
the table and the figure cannot disagree about what was dropped.

### There is no residual airmass dependence

Five within-setup visit pairs, where central wavelength, star and slicer are all
held fixed and airmass is the only difference:


| setup      | @6400    | @6600    |
| ---------- | -------- | -------- |
| 2023-10-15 | -2.4 pts | -1.9 pts |
| 2023-10-16 | **+1.5** | **+2.3** |
| 2023-10-17 | **+1.1** | **+1.6** |
| 2023-10-20 | -1.6     | -2.5     |
| 2023-11-08 | —        | -0.9     |


Signs disagree, mean about -0.3 pts. **The extinction correction works.**
An earlier claim of a -4.5%/airmass trend, made on one pair, is withdrawn.

### The scatter is photons, not processing

Comparing the star's counts per second — extinction removed, so visits at
different airmass are comparable — against the throughput derived from them:

```
correlation(counts, throughput) = +0.98   (n=8)
spread in counts      5.82%
spread in throughput  5.20%
```

**Read this carefully, because it proves less than it first appears to.** All 8
points share one star and one central wavelength, so the expected photon count is
a fixed number and throughput is proportional to counts *by construction*. The
correlation could hardly come out otherwise.

What the numbers do establish is narrower:

- The correlation is **+0.98 rather than +1.00**, and the two spreads nearly
match (5.82% against 5.20%). So the sensfunc fit, the trim interval and the
zeropoint polynomial add only a little independent wobble — the scatter is
**inherited from the extracted counts**, not manufactured downstream of them.
- Low counts at fixed airmass and fixed configuration means either fewer photons
arrived or the extraction failed to collect them. The enclosed-flux measurement
bounds the second: at a common 3.4" aperture these 8 cubes enclose 0.94-0.96,
except 2023-10-20 B v2 at 0.81. So aperture losses cannot account for a 5.8%
spread.

Sharing one central wavelength, one star and one slicer also excludes blaze
mixing and the choice of standard. Star concentration correlates at +0.00.

What remains is transparency, telescope or instrument state — **not separable
from this dataset**. Within-night variation (1.47 pts) is as large as
night-to-night (1.21 pts). 2023-11-07 shows the effect directly: three
consecutive 60 s frames varying 18% in counts at constant airmass, with the star
broadening as counts fell.

### Comparing RH1 with RL is not one number against one number

Each curve peaks near its own central wavelength — the blaze tracking the grating
angle. RH1 spans six central wavelengths from 6140 to 6680, so the product is
**throughput as a function of wavelength with the central wavelength named**, not
a grating-wide average. A value quoted at a fixed wavelength across different
central wavelengths mixes blaze positions.

The two former low outliers, **2024-03-13** and **2023-11-07**, are no longer
treated as outliers to be explained away — they are excluded, on the config and
achromatic-ratio evidence above.

---



## 9. Tools


| file                                              | purpose                                                                               |
| ------------------------------------------------- | ------------------------------------------------------------------------------------- |
| `pypeit_test/build_rh1_template.py`               | stitch or cut a template; `--only`/`--max-seeds 1`/`--target` for one-exposure builds |
| `pypeit_test/gate_rh1_template.py`                | prove a template on a real setup before trusting it                                   |
| `pypeit_test/pick_seed_slits.py`                  | rank slices; separate seed quality from setup health                                  |
| `pypeit_test/check_slice_wavelength_agreement.py` | find slices shifted from their neighbours                                             |
| `pypeit_test/check_seed_agreement.py`             | cross-correlate two independent solutions where they overlap                          |
| `pypeit_test/check_skyline_wavelengths.py`        | verify the scale against sky lines the arc never saw                                  |
| `pypeit_test/plot_rh1_throughput.py`              | all throughput curves on one axis, RL underneath                                      |
| `pypeit_test/plot_slice_diagnostics.py`           | which slice is wrong, wavelength and geometry                                         |
| `pypeit_test/plot_1108_slice_defect.py`           | the dead column, from the raw flat                                                    |
| `pypeit_test/rh1_throughput_table.py`             | throughput table, edge-flagged, with visit-pair comparison; holds `DROP`              |
| `RH1 pypeit run/run_rh1.py`                       | night-by-night driver to spec2d                                                       |
| `RH1 pypeit run/run_rh1_throughput.py`            | cubes to throughput, split by airmass; `--boxcar`/`--skip-coadd`                      |
| `RH1 pypeit run/redo_6520_setups.py`              | re-reduce a setup against a narrower template                                         |


Templates live in `pypeit_test/` rather than the conda environment, so an
environment rebuild cannot silently lose them. PypeIt accepts an absolute path
for `reid_arxiv`.

**Calibrations are keyed by frame, not by parameters.** A `Calibrations/`
directory built before a parameter change is *reused*, not rebuilt. Delete it
whenever the template or a slit-edge parameter changes, or the old solution
survives into the science run silently.

---



## 10. Still open

- **Two cubes lose flux that the aperture fix did not recover.** At a common
3.4" aperture, **2023-10-20 B v2** encloses 0.81 and **2023-11-08 A v2** 0.84,
where their own siblings — same night, same setup, same calibrations — reach
0.94-0.96. Not weather, since a sibling visit hours away is fine; not aperture,
since the aperture is now identical. Unexplained. Both flag in the table.
- **Column 653** may be dead on other KCRM data. If it is a permanent detector
defect it belongs in the bad-pixel mask for all reductions, not patched per
night. One night is not enough to tell.
- **Separating transparency from telescope state** needs an external handle —
for instance, whether RL and RH1 measurements sharing a date move together.
The two excluded nights are the obvious test case.
- Whether the two standards agree: with 2024-03-13 excluded, **feige34
contributes two usable points** (30.6% at 6200 A, 32.4% at 6600 A) against
g191b2b's 31.7% (n=1) and 29.3% (n=10). The 6200 A pair agrees to 1.1 pts;
the 6600 A pair differs by 3.1 pts, with feige34 above every one of the ten
g191b2b values. Suggestive of a star-dependent offset, not yet conclusive on
two points.
- **Only 2 cubes now carry @6200 A.** The number is tight (+-0.56%) precisely
because the non-photometric nights left, but two points constrain a mean, not
a scatter. Any future RH1 blue-end data should go here first.

