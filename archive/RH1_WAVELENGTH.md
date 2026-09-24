# RH1 needs a wavelength template PypeIt does not ship

Why the RH1 reduction cannot start the way the RL one did, and how the missing
piece is built. Written 2026-08-24, against PypeIt 2.0.1.

Scope: 9 nights, 10 configs, 44 standard-star frames. The RH1 standards are
**g191b2b** (36 frames, 7 nights) and **feige34** (8 frames, 2 nights) — not the
feige110 of the RL run. Both are calspec standards and both sit within 2" of
their catalogue positions, so `pypeit_sensfunc` will match them on coordinates
as it did for feige110. Nothing in this document depends on which star it is:
the problem below is a property of the grating.

---

## Vocabulary

| term | meaning |
|---|---|
| **arc frame** | exposure of the FeAr calibration lamps; a comb of emission lines at known wavelengths |
| **line list / catalogue** | table of true wavelengths for those lamp lines (here `FeI`, `ArI`, `ArII`) |
| **wavelength solution** | the polynomial turning detector pixel into wavelength, for one slice |
| **reid_arxiv / template** | a stored reference arc spectrum *with* its wavelength solution, shipped with PypeIt per grating |
| **`full_template`** | method that matches an observed arc against that template to decide which catalogue line each feature is |
| **holy-grail** | fallback method that tries to identify lines from scratch, with no reference |
| **slice** | one of the 24 strips the image slicer cuts the field into; each gets its own solution |
| **central wavelength** | the wavelength at the middle of the detector, set by the grating angle (header `RCWAVE`) |

---

## The failure

The first night, 2023-09-16, died in calibration:

```
ValueError: zero-size array to reduction operation minimum which has no identity
  flatfield.py:2220  mnmx_wv[slit_idx, 0] = np.min(waveimg[onslit_init])
```

Read backwards, that is a slice with no valid wavelengths at all. The cause is
one line further up the log:

```
keck_kcwi.py:config_specific_par  Full template solution is unavailable
                                  Adopting holy-grail algorithm
```

PypeIt assigns a template per grating, and its list has a hole:

| grating | template |
|---|---|
| BL, BM, BH2, BH3 | `keck_kcwi_*.fits` |
| RL, RM1, RM2, RH3 | `keck_kcrm_*.fits` |
| **RH1, RH2, RH4** | **none — falls back to holy-grail** |

Not a local problem: no `keck_kcrm_RH1.fits` exists in PypeIt's `reid_arxiv`
directory on the `develop` branch either. RL worked precisely because it is on
the supported list.

Holy-grail then solved **1 slice of 24** on this setup:

```
slit  6 (spat  515): rms=  0.346 px  nlines= 59  5863-6495 A     <- usable
slit  7 (spat  588): rms=  0.129 px  nlines=  6  -84441-23817 A  <- "fit", nonsense
slit  9 (spat  736): rms=  0.192 px  nlines=  6  -51251-113865 A <- "fit", nonsense
the other 21                                     no fit
```

Note slits 7 and 9: a low RMS through six points spanning 10^5 A. A returned fit
object is not a solution, which is why `pick_seed_slits.py` also demands line
count, physical span and smooth monotonic dispersion.

The 21 unsolved slices are what the flat field then trips over. Every RH1 night
would fail identically, so the driver was stopped after one night.

---

## Why one setup cannot supply the template

RH1 is a high-dispersion grating: one exposure covers only about **630 A**. The
nine selected nights were taken at six different central wavelengths, so their
arcs sit at different places on the spectrum:

```
      5800      6000      6200      6400      6600      6800      7000 A
        |         |         |         |         |         |         |
6140  [=========================]                              2023-11-17
6150   [=========================]                             2024-03-13
6200     [=========================]                           2023-09-16, 11-07
6520              [=========================]                  2023-10-15/16/17/20
6600                 [=========================]               2023-11-17 (2nd cfg)
6680                    [=========================]            2023-11-08 (1x1)
        |<------------------ 5824-6996 A needed ------------->|
```

A template must cover the union, so it is **stitched** from the best slice of
several setups — which is what `pypeit.core.wavecal.templates.build_template`
is for. Consecutive seeds are cut at the midpoint of their overlap, and
`shift_wave=True` aligns each new snippet to the one before it.

---

## How much must a template's own scale be trusted?

Less than it first appears, and this is the load-bearing point.

`full_template` uses the template **only to decide which catalogue line each
observed feature is**. The wavelength solution it then reports is fit to
*catalogue* wavelengths, not to the template. So an error in the template's own
scale does not propagate into the science; it only has to be small enough to
name each line correctly — a fraction of the **line spacing**, here about
4.2 A, with PypeIt's `match_toler` at 2.0 A.

That matters because the two seeds available first do *not* agree perfectly.
Cross-correlating their arc spectra where they overlap (266 A, independently
solved, different nights and different slices):

```
cross-correlation peak 0.895 at +1.131 A   (0.27 line spacings)
zero-shift correlation 0.190
```

A real 1.13 A systematic, in the region where each solution is weakest — it is
the red end of one and the blue end of the other, where each has the fewest
fitted lines. It is under the 2.0 A that would cause a misidentification, but
the margin is thin, and it is *not* visible in either fit's own RMS (0.35 and
0.50 px): a fit only scores itself against the lines it chose.

`check_seed_agreement.py` performs this test. It deliberately does not compare
the two solutions' `wave_fit` arrays — those are catalogue values, identical
whenever both picked the same line, and blind to a solution that put that line
at the wrong pixel.

---

## The gate

The seed comparison is a warning light, not the verdict. The verdict is
empirical: build the template, run `full_template` with it on the night that
crashed, and count how many of the 24 slices solve and at what RMS.
`gate_rh1_template.py` does exactly that, on 2023-09-16 setup B:

| | slices solved | median RMS | max RMS | run_pypeit |
|---|---|---|---|---|
| holy-grail | **1 / 24** | 0.346 px | — | crashed in the flat field |
| `full_template` | **24 / 24** | **0.123 px** | 0.244 px | returned 0 |

Every slice was fit with 55-69 catalogue lines, and the fitted dispersion agrees
across all 24 to 0.3% (0.3050-0.3059 A/px) — an internal consistency check the
holy-grail run could not offer at all, having produced one solution.

So the 1.13 A seed disagreement did not matter, exactly as the argument above
predicts: it is under half a line spacing, the lines were named correctly, and
the reported solution is a fresh fit to catalogue wavelengths. The template's
job is identification, and it did it.

Note the template is now *better than the seeds it was built from* — 0.123 px
median against the seed's 0.346. That is not circular: correct identification is
all the template contributes, and the fit is to the catalogue either way. It
does mean the template should be rebuilt from a passing night's solutions, which
`build_rh1_template.py` does automatically by ranking every slice of every
WaveCalib it can find.

### Growing the coverage

The 24 slices of one setup do not share a wavelength window — the slicer feeds
the grating at a slightly different angle per slice, spreading the windows over
about 40 A (5860.1 to 5899.3 at the blue end here). So the union across slices
reaches slightly further than any single slice, and the builder picks the bluest
and reddest rather than only the lowest-RMS.

Beyond that, coverage grows by bootstrapping: a template covering most of a new
setup's range solves that setup, and its solutions then extend the template.
2023-11-17 is the useful one, carrying both a 6140 and a 6600 config — the
bluest and nearly the reddest of the nine nights — each about 90% covered by the
first template. Both solved completely on that partial coverage:

| setup | central λ | slices | median RMS | span solved |
|---|---|---|---|---|
| 2023-09-16 B | 6200 | 24/24 | 0.123 px | 5860.1-6526.9 A |
| 2023-10-15 B | 6520 | 24/24 | 0.064 px (best slice) | 6194.5-6823.9 A |
| 2023-11-17 B | 6140 | 24/24 | 0.114 px | 5796.5-6471.0 A |
| 2023-11-17 C | 6600 | 24/24 | 0.135 px | 6278.3-6899.2 A |

2023-10-15 is the clearest demonstration: holy-grail had given it 1 slice of 24
at RMS 0.499, and rerunning the identical frames against the template gave 24 of
24 with a best slice at 0.064 — the same data, only the line identification
changed.

The template in use is rebuilt from the two configs that reach furthest in each
direction, both now from validated 24/24 solutions rather than holy-grail seeds:

```
5796.5 A |======= 2023-11-17 B (6140) =======|                    6431.7
                              6311.9 |======= 2023-11-17 C (6600) =======| 6899.2
                                  stitched at 6371.8 (120 A overlap)
```

3720 pixels, 5796.5-6899.1 A, `BINSPEC = 2`, monotonic. The dispersion runs
0.2878 A/px at the blue end to 0.3056 at the red — the grating's own trend, not
a step at the joint.

That covers every selected setup except the red end of 2023-11-08 (6680, which
needs about 6975 A), leaving ~76 A — 12% of that setup's range — to be
extrapolated. Both 2023-11-17 configs solved on a comparable shortfall, so this
is expected to hold, and the driver reports the slice count per setup either way.

2023-11-08 is also the only 1x1 setup, so its template is resampled by
`full_template` — a path PypeIt's own comment marks "not yet tested". Checked
directly: 3233 -> 6466 pixels, wavelength endpoints identical to 4 decimals,
dispersion exactly halved, flux peak preserved. That path is sound; the setup
fails for a different reason.

---

## Why 2023-11-08 needs its own template

It is the only **Small** slicer setup. Slicer width sets spectral resolution, so
its arc lines are physically narrower than anything the template was built from,
and a template cannot identify lines it does not resolve.

**Measure line width in Angstroms, never in pixels.** In pixels this setup and a
known-good RL one look identical — both ~2.7 px — which is what made this take
three attempts to pin down. The pixel is a different physical size in each case:

| spectrum | FWHM px | A/px | **FWHM (A)** | template/data |
|---|---|---|---|---|
| shipped `keck_kcrm_RL` template | 4.00 | 1.8356 | 7.34 | — |
| RL 2023-12-17, Small, 1x1 — **works** | 6.50 | 0.9181 | 5.97 | **1.23x** |
| our RH1 template (Large seeds) | 5.00 | 0.2973 | 1.42 | — |
| RH1 2023-11-17 C, Large, 2x2 — works | 5.00 | 0.2841 | 1.42 | 1.00x |
| RH1 2023-11-08, Small, 1x1 — **fails** | 4.50 | 0.1531 | **0.69** | **2.06x** |

The RL Small-slicer night works because RL is low-dispersion: its lines are
5.97 A wide against a 7.34 A template, near-matched. RH1 at the Small slicer
reaches 0.69 A against our 1.42 A template, so the data resolves blends the
template holds as single features — the line inventories differ.

### What this is not

Two explanations were tested and ruled out by measurement, both worth recording
because each looked convincing:

**Not wavelength coverage.** The uncovered red end is 12%, and the fully covered
blue half failed too. 2023-10-15 was 100% covered, scored a poor global cc of
0.34, and still solved 24/24.

**Not template alignment.** A template much wider than one exposure does make
alignment ambiguous — `full_template` cuts a window as wide as your data out of
the template at a position chosen by one cross-correlation:

```python
i0    = npad // 2 + int(shift_cc)
mspec = temp_spec[i0 : i0 + nspec]
```

Rebuilding at one exposure wide (`--only <setup> --max-seeds 1 --target LO HI`)
cut the shift from 2065 px to 261-493 px and raised the local correlation from
cc 0.34 to a median of 0.767 — and the setup still solved 0 of 25. Alignment was
real but never the blocker.

The fix is therefore a separate Small-slicer template, seeded by holy-grail from
this setup's own arcs, exactly as the Large-slicer template was originally
seeded.

---

## Checking the science, not just the calibrations

Wavelength solutions being sound does not make the reduction sound, so the first
finished frame was checked with `check_skysub.py` (2023-09-16, g191b2b):

```
STD_CHIS          0.95 - 1.03 across the 24 slices   (1.0 is a correct noise model)
star located on   slices 11-13, peaking at slice 12  (matches user_regions :40,66:)
verdict           SKY MODEL SUSPECT on all 24 slices
```

The verdict is a false alarm, and the control proves it: the same tool on
2024-12-04 B — the RL night behind the published 29.5% throughput — returns
`SKY MODEL SUSPECT` on most of its slices too, at STD_CHIS 0.955-0.976 against
RH1's 0.95-1.03. The verdict fires when the sky-zone residual exceeds 5% of the
sky level, which on a short standard-star exposure it routinely does simply
because the sky level is low. **On standard-star frames read STD_CHIS, not the
verdict.** For scale, the RL skysub investigation measured STD_CHIS 3.235 for a
genuinely corrupted sky model and 1.322 for the accepted `joint_fit` fix.

## Tools

| file | purpose |
|---|---|
| `pypeit_test/pick_seed_slits.py` | rank a WaveCalib's slices; reject fits that are nonsense despite low RMS |
| `pypeit_test/check_seed_agreement.py` | cross-correlate two independent solutions where they overlap |
| `pypeit_test/build_rh1_template.py` | stitch the seeds into `keck_kcrm_RH1.fits` |
| `RH1 pypeit run/run_rh1.py` | night-by-night driver; injects the template and the sky regions |

The driver writes into each setup's `.pypeit` file:

```ini
[calibrations]
    [[wavelengths]]
        method = full_template
        reid_arxiv = <repo>/pypeit_test/keck_kcrm_RH1.fits
```

An absolute path is accepted: PypeIt's `get_file_path` returns any path that
already resolves, so the template stays in the repo rather than in the conda
environment, where an environment rebuild would silently lose it.

**Calibrations are keyed by frame, not by parameters.** A `Calibrations/`
directory built before the template existed would be *reused*, not rebuilt, and
the holy-grail solution would survive into the science run. The driver deletes
any such directory when it inserts the template.
