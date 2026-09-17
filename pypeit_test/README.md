# PypeIt test case — feige110, night 2023-12-17, KCWI red

A sandbox for working through the PypeIt tutorial
(https://pypeit.readthedocs.io/en/latest/tutorials/kast_howto.html) on real data
from this project. Nothing here touches the raw archive: `raw/` holds **symlinks**
into `../fits/`, and PypeIt only ever writes into `keck_kcrm_*/`.

## What is in raw/

One observing night, chosen because every calibration comes from that same night:

    4  science    feige110, RL grating, Small IFU, 1x1 binning
    6  arclamp    ThAr + FeAr
    6  flatlamp   continuum
    6  domeflat
    3  contbars
    3  twiflat
    7  bias       1x1, AMPMODE L2U2 — fetched separately, see below
    1  bias       2x2 FPCam TEST — junk, PypeIt isolates it as Setup A

Science and dispersed calibrations come from `../outputKC/same_night_pairs.csv`,
night `2023-12-17`. The seven usable biases do not, and could not: see below.

## Running it

    conda activate pypeit
    cd /Users/kaiwenzhang/PycharmProjects/KCWI/pyKOA/pypeit_test

    pypeit_setup -r raw -s keck_kcrm -c all -d .    # writes keck_kcrm_A and _B
    cd keck_kcrm_B
    run_pypeit keck_kcrm_B.pypeit                   # ~2h50m for these 4 frames

Use `keck_kcrm_B`, not `A` — see below. `run_pypeit -o` overwrites outputs in
place. To start over from nothing (pypeit_setup will not overwrite existing
directories):

    rm -rf keck_kcrm_A keck_kcrm_B *.log

## Results of the first run

Completed, exit 0, in 2h50m. Four `spec2d` files, ten calibration frames, QA in
`QA/MF_B.html` — an index of the 302 PNGs in `QA/PNGs/`; open it in a browser
(`open QA/MF_B.html`), not in the editor.

No `spec1d`, which is expected for an IFU. Extraction happens after building a
cube:

    pypeit_coadd_datacube      # spec2d -> datacube
    pypeit_extract_datacube    # datacube -> 1D spectrum

Non-fatal warnings worth checking in QA: ~100 poorly-conditioned polynomial fits
and 57 `maxiter` hits in the wavelength/tilt solutions (`QA/PNGs/Arc_1dfit_*`),
slit-edge shifts of 0.6-1.0% indicating mild within-night flexure, and no tuned
scattered-light model for grating RL (it fitted from generic defaults and
converged).

## Three things to know

**Use spectrograph `keck_kcrm`, not `keck_kcwi`.** KCRM is the red channel; these
are `KR.*` files. Setup prints a warning that the headers say `INSTRUME = KCWI`
while KCRM was expected -- that is how KOA labels the instrument as a whole, and
the frame typing is unaffected.

**`-c A` is the wrong configuration.** PypeIt found two, and sorts them by its own
rules, so the useful data is not first:

    Setup A:  binning 2,2   decker FPCam  dispname None    1 file
    Setup B:  binning 1,1   decker Small  dispname RL     35 files

**The biases cannot come from the caliblist, by design.** KCWI tags a properly
taken bias sequence `stateid = 0`, `statenam = "Shared Bias"` -- correctly, since a
bias depends only on detector readout and is shared across every instrument
configuration. But KOA's `nph-getCaliblist` associates calibrations *by stateid*,
so a `stateid = 0` bias matches no science frame and is never returned. 68% of all
red biases in KOA (3517 of 5144) are labelled this way.

What the caliblist does return is the accident: an engineering `TEST` readout
taken at 00:10 while the night's state was still active inherits that stateid and
matches -- despite being 2x2 with the decker parked at FPCam, useless for 1x1
science. That is the frame PypeIt quarantines as Setup A.

The seven real biases were found by querying the archive directly, keyed on the
two UT dates the Keck night spans and filtered to the science frame's readout:

    select koaid, instrume, filehand, binning, ampmode from koa_kcwi
    where koaimtyp='bias' and camera='RED' and binning='1,1'
      and koaid like 'KR.20231218.08%'

They match the science frames on BINNING, AMPMODE, CCDMODE, GAIN1, CDSSPEED and
array shape (4234x4188). Note `date_obs` is an Oracle DATE column, so string
comparison on it fails with ORA-01861; the date embedded in `koaid` works.

The three twilight flats are commented out in the .pypeit file: PypeIt could not
assign them a frametype. That is expected and harmless -- uncomment and label them
`illumflat` only if you want them in the illumination correction.
