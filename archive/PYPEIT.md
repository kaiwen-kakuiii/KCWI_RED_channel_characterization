# Reducing a night with PypeIt

Picks up where [download/DOWNLOAD.md](../download/DOWNLOAD.md) stops: you have a night in
`fits/by_night/<night>/` and want a flux-calibrated spectrum out of it.

Written against **PypeIt 2.0.1**, spectrograph **`keck_kcrm`** — KCRM is the red
channel, and these are `KR.*` files. Setup will warn that the headers say
`INSTRUME = KCWI` while KCRM was expected; that is how KOA labels the instrument
as a whole, and frame typing is unaffected.

Worked example throughout: **feige110, night 2023-12-17** — 4 science frames plus
31 calibrations, the night the pipeline's own tests use.

```
fits/by_night/2023-12-17/          35 files
        |  pypeit_setup
        v
keck_kcrm_A.pypeit                 one setup, frametypes assigned
        |  run_pypeit                          ~2h50m for 4 frames
        v
Science/spec2d_*.fits              4 files, ~620 MB each
Calibrations/*.fits                10 products
        |  pypeit_coadd_datacube
        v
feige110_2023-12-17.fits           the datacube
        |  pypeit_extract_datacube
        v
spec1d_*.fits
        |  pypeit_sensfunc
        v
feige110_sens.fits                 the sensitivity function
```

Each stage below has a **Check** section. Run them. A wavelength solution that
converged is not the same as a wavelength solution that is right.

---

## 1. Setup

```bash
pypeit_setup -s keck_kcrm -r fits/by_night/2023-12-17 -c all -d pypeit_run
cd pypeit_run/keck_kcrm_A
```

One night per `-r`. Pointing it at several nights at once is the mistake the
per-night tree exists to prevent: PypeIt groups by `dispname`, `decker`,
`binning`, `cenwave` and no part of that key is a date, so it will happily build
one master arc out of four nights of arcs.

**Check the grating has a wavelength template first.** PypeIt ships one for RL,
RM1, RM2 and RH3 only. For RH1, RH2 and RH4 it prints `Full template solution is
unavailable`, falls back to holy-grail, and — measured on RH1 — solves 1 slice of
24, after which the flat field dies on the other 23.

If you are reducing a grating other than RL, read
**[RH1_PROCEDURE.md](RH1_PROCEDURE.md)** before starting. It is the delta from
this document: the missing template and how to build one, why a second template
is needed per slicer, why a template must be one exposure wide, a dead detector
column that invents a slice, why the extraction aperture must be fixed in arcsec
rather than left to its default, and the extra cube-splitting step.
[RH1_WAVELENGTH.md](RH1_WAVELENGTH.md) carries the template detail.

### Check

Read the file before running anything.

```bash
grep -A6 "setup read" keck_kcrm_A.pypeit
awk -F'|' '/\.fits/{gsub(/ /,"",$2); print $2}' keck_kcrm_A.pypeit | sort | uniq -c | sort -rn
```

**One** setup, 35 files, and this typing:

| count | frametype | from |
|---|---|---|
| 7 | `bias` | the night's shared bias set |
| 6 | `pixelflat,scattlight` | flatlamp |
| 6 | `illumflat,trace` | domeflat |
| 4 | `science` | feige110 |
| 3 | `arc` + 3 `tilt` | the 6 arclamps, split |
| 3 | `align` | contbars |
| 3 | `None` | twiflat — commented out |


The three twilight flats stay commented: PypeIt cannot assign them a frametype.
Harmless. Uncomment and label them `illumflat` only if you want them in the
illumination correction.

---

## 2. Reduce

```bash
run_pypeit keck_kcrm_A.pypeit
```

About 2h50m for these four frames. `run_pypeit -o` overwrites outputs in place.
`pypeit_setup` will not overwrite an existing directory, so to start clean:
`rm -rf pypeit_run`.

### Check — did it finish?

```bash
tail -3 keck_kcrm_A.log                                   # "Data reduction complete"
grep -E "^\s*(Traceback|.*\bERROR\b)" run.log | grep -v "centroid error"
```

The `grep -v` is not optional. `Max centroid error: None` is an INFO line that a
naive search for ERROR reports as a failure — it is the only thing a clean run
here matches.

### Check — the calibrations

```bash
pypeit_chk_edges      Calibrations/Edges_A_1_DET01.fits.gz
pypeit_chk_flats      Calibrations/Flat_A_1_DET01.fits
pypeit_chk_wavecalib  Calibrations/WaveCalib_A_1_DET01.fits
pypeit_chk_tilts      Calibrations/Tilts_A_1_DET01.fits
pypeit_chk_scattlight Calibrations/ScatteredLight_A_1_DET01.fits.gz
pypeit_chk_alignments Calibrations/Alignment_A_1_DET01.fits
```

`chk_edges` and `chk_tilts` drive **ginga**; the rest use matplotlib.

`chk_wavecalib` is the one to weigh most — it prints per-slit RMS. Measured on
this night: ~100 poorly-conditioned polynomial fits and 57 `maxiter` hits in the
wavelength and tilt solutions. They are warnings, not failures, but "it
converged" is not evidence the solution is good. Look at `QA/PNGs/Arc_1dfit_*`
and judge the RMS.

Two other known warnings for this night, both non-fatal: slit-edge shifts of
0.6–1.0%, which is mild within-night flexure, and no tuned scattered-light model
for grating RL — it fits from generic defaults and converges.

### Check — the 2D spectra

```bash
pypeit_show_2dspec      Science/spec2d_KR.20231218.17261.35-*.fits
pypeit_chk_noise_2dspec Science/spec2d_*.fits
pypeit_parse_slits      Science/spec2d_KR.20231218.17261.35-*.fits
```

`chk_noise_2dspec` histograms (data − model)/σ. It should look like a unit
Gaussian. A wider one means the noise model is wrong, and so is every error bar
downstream of it.

QA index: `open QA/MF_A.html` — in a browser, not the editor. It indexes ~300
PNGs.

There is no `spec1d` at this stage and that is correct for an IFU. Extraction
happens after the cube.

---

## 3. Datacube

Write a `.coadd3d` file next to the spec2d files:

```ini
[rdx]
    spectrograph = keck_kcrm

[reduce]
    [[cube]]
        combine = True
        output_filename = feige110_2023-12-17.fits
        save_whitelight = True

spec2d read
filename
Science/spec2d_KR.20231218.17261.35-feige110_KCRM_20231218T044741.338.fits
Science/spec2d_KR.20231218.17512.28-feige110_KCRM_20231218T045152.243.fits
Science/spec2d_KR.20231218.17770.37-feige110_KCRM_20231218T045610.406.fits
Science/spec2d_KR.20231218.18023.38-feige110_KCRM_20231218T050023.386.fits
spec2d end
```

```bash
pypeit_coadd_datacube feige110.coadd3d
```

`combine = True` is right for this night because all four exposures share one
pointing — identical RA/DEC, `ra_off = dec_off = 0`. There is no dither, so no
alignment is needed. Check that before copying this for another night: if the
frames are offset, set `align = True`.

Validate the file before spending time on it:

```bash
python -c "from pypeit.inputfiles import Coadd3DFile; Coadd3DFile.from_file('feige110.coadd3d')"
```

It prints `.coadd3d file successfully vetted.` and resolves every path.

### Check

```bash
pypeit_view_fits keck_kcrm feige110_2023-12-17_whitelight.fits
```

That is what `save_whitelight` is for. The star should be a single compact
source near the field centre. Smeared or doubled means the exposures did not
stack — set `align = True` and redo.

Then confirm the cube's shape:

```bash
python -c "
from astropy.io import fits
h = fits.open('feige110_2023-12-17.fits'); h.info()"
```

---

## 4. Extract the 1D spectrum

```bash
pypeit_extract_datacube feige110_2023-12-17.fits
```

`-b BOXCAR_RADIUS` widens the aperture. For a standard you want essentially all
the flux, so err large.

### Check

```bash
pypeit_show_1dspec      spec1d_*.fits
pypeit_chk_noise_1dspec spec1d_*.fits
```

A smooth continuum across the grating's range with no steps at slice boundaries,
and Balmer absorption where feige110 — an sdO — should show it. Steps at slice
boundaries point back at the illumination correction, which is where the
twilight flats would come in.

---

## 5. Sensitivity function

```bash
pypeit_sensfunc spec1d_*.fits -o feige110_sens.fits --algorithm UVIS
```

No manual star naming. PypeIt matches on coordinates, and the science header
lands **0.30″** from calspec's `feige110_stisnic_008.fits.gz` entry
(23:19:58.4 −05:09:56.17), well inside any tolerance.

`--algorithm UVIS` suits KCRM at 8150 Å centre; `IR` is for the near-infrared
where telluric fitting dominates. Check which one the default picks before
accepting the result.

### Check

`pypeit_sensfunc` writes its own QA. The zeropoint should vary smoothly with
wavelength. Sharp features usually mean the extraction pulled in a cosmic ray,
or the standard's own absorption lines were divided out badly.

Comparing this zeropoint against other nights and gratings is the red-channel
characterisation this repo is built for.

---

## Tool reference

Every check tool in PypeIt 2.0.1, and what it takes:

| Tool                      | Argument                                | Display    |
| ------------------------- | --------------------------------------- | ---------- |
| `pypeit_chk_edges`        | `Edges_*.fits.gz`                       | ginga      |
| `pypeit_chk_flats`        | `Flat_*.fits`                           | matplotlib |
| `pypeit_chk_wavecalib`    | `WaveCalib_*.fits` or a spec2d          | matplotlib |
| `pypeit_chk_tilts`        | `Tilts_*.fits`                          | ginga      |
| `pypeit_chk_scattlight`   | `ScatteredLight_*.fits.gz`              | matplotlib |
| `pypeit_chk_alignments`   | `Alignment_*.fits`                      | matplotlib |
| `pypeit_chk_flexure`      | spec2d or spec1d                        | matplotlib |
| `pypeit_show_2dspec`      | a spec2d file                           | ginga      |
| `pypeit_show_1dspec`      | a spec1d file (not coadd_1dspec output) | ginga      |
| `pypeit_chk_noise_2dspec` | spec2d file(s)                          | matplotlib |
| `pypeit_chk_noise_1dspec` | spec1d file(s)                          | matplotlib |
| `pypeit_parse_slits`      | a spec2d or `Slits_*` file              | text       |
| `pypeit_view_fits`        | `<spectrograph> <file>`                 | ginga      |
| `pypeit_qa_html`          | rebuilds the QA index                   | —          |

---

## Gotchas

**The setup letter is not stable.** PypeIt names configurations in its own order,
so the useful data is not always `A`. An earlier hand-built run of this night
included the 2×2 FPCam test bias, which PypeIt isolated as Setup A and pushed the
real 35 files into Setup B. The pipeline no longer downloads that frame, so there
is now one setup — but always read the setup block rather than assuming a letter.

**`pypeit_setup` will not overwrite.** It refuses to write into an existing
directory. `rm -rf pypeit_run` to start over.

**QA HTML can be incomplete.** One run here wrote 3 of 4 per-frame QA files, the
log ending mid-generation. The science products were all present and correct.
`pypeit_qa_html` rebuilds the index if you want it whole.

**No `spec1d` after `run_pypeit`.** Expected for an IFU. If you are looking for
one before building the cube, you are looking too early.
