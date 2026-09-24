# KCRM red-channel throughput

Measured throughput of the seven KCRM gratings, from KOA standard-star frames
reduced with PypeIt. Raw frames and PypeIt products are not in this repo. They
are fetched and reduced locally, and the published curves are what is kept.

**Download.** [download/DOWNLOAD.md](download/DOWNLOAD.md) finds the
standard-star frames and writes `fits/by_night/<grating>/<night>/`.

**Reduce a new night.** [pipeline/TUTORIAL.md](pipeline/TUTORIAL.md). The
grating files are in `pipeline/config/`. Wavelength templates and the ThAr
line lists are in `templates/`.

**Results.** [results/GRATING_SUMMARY.md](results/GRATING_SUMMARY.md) is the
register. Per-grating curves are `results/<GRATING>/<night>/keck_kcrm_<letter>/`.
Each night keeps the published `*_sens_IR_p15cov.fits` or `*_sens_IR_p5cov.fits`
file, the `.sens` recipe, the `.coadd3d` frame list, and the `.pypeit` file.
Figures and tables sit next to those folders. Rewrite them with:

```bash
python pipeline/utils/plot_blaze_summary.py
python pipeline/utils/blaze_tables.py --render
```

`blaze_tables.py` reads `results/grating_summary.template.md`. Do not edit
`results/GRATING_SUMMARY.md` by hand.

**Archive.** [archive/](archive/) holds the writeups and the drivers that
produced the published nights.
