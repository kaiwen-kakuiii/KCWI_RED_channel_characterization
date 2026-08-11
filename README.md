# KCWI red-channel standard stars

Three scripts that find the KCWI/KCRM standard-star frames in the Keck archive,
tell you what coverage they give, and download them together with calibrations
taken on the same night.

```
kcwi_stds_selected.txt          40 standards, an observer's starlist
        |
        |  fetch_std_red.py     one ADQL cone query, label each frame with its star
        v
outputKC/std_red_matched.csv    1307 red science frames, 27 standards
        |
        |------> analyze_gratings.py    read-only: which star covers which grating
        |
        |  download_std_red.py  filter, gate on same-night calibs, download
        v
fits/lev0/                      science frames
fits/calib_same_night/lev0/     their calibrations
fits/by_night/<night>/          the same files, one directory per night  -> PypeIt
```

Run `fetch_std_red.py` from the repo root, the other two from `pyKOA/`. Paths
inside `download_std_red.py` resolve against the script, so it works from either.

---



## 1. `fetch_std_red.py` — find and label

Asks the archive once for every red science frame that falls inside any
standard's cone, then works out locally which standard each frame belongs to.

No command-line options. Three constants at the top:


| Constant       | Default                  | What it does                                                                                                                                                               |
| -------------- | ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `STARLIST`     | `kcwi_stds_selected.txt` | the standards to look for                                                                                                                                                  |
| `MATCH_RADIUS` | `0.05 deg`               | server-side cone. Must stay well wider than the IFU: it has to catch frames taken at a star's proper-motion-corrected position, up to ~30″ from its epoch-2000 coordinates |
| `FOV_RADIUS`   | `19.4″`                  | local test for whether the star was actually *in* the IFU. Half-diagonal of the largest slicer, 33″ × 20.4″                                                                |


**How it decides.** One query ORs all 40 cones into a single `WHERE`, so the
server does the spatial cut — one request, ~9 s, ~600 KB back. Each frame is
then labelled with its nearest standard and kept only if both:

- `in_fov` — the star was inside the 19.4″ field, testing against **both** the
epoch-2000 and the proper-motion-corrected position and accepting whichever is
closer, because observers were inconsistent about correcting at the telescope
- `name_ok` — the observer's `TARGNAME` names that standard, after normalising
(`bd02d3375`, `BD023375` and `BD+02 3375` all collapse to the same string)

**What it deliberately does not do:** filter on the `object` string. `OBJECT` and
`TARGNAME` are sticky — whichever the observer last set persists into following
exposures — so a frame labelled `DOME/Dark/240.0 sec.` is often a real on-target
standard exposure carrying a stale label. 29 of 40 such frames have an exposure
time contradicting their own label. The cone cut is the trustworthy selector.

Re-running is cheap: the raw query is cached as
`kcwi_red_cones_<starlist>.tbl` and only re-fetched if you delete it. The cache
is named after the starlist because the labelling is only valid while the cones
that produced it match the stars being labelled against.

---



## 2. `analyze_gratings.py` — inspect coverage

Read-only. Prints to stdout, writes nothing, downloads nothing. Run it to decide
which star and grating to actually work on.

```
python analyze_gratings.py
```

For all IFUs combined and then for Large, Medium and Small separately:

- stars ranked by how many distinct gratings they were observed with, with the
frame count on each
- total frames per grating
- the distinct `ifunam` values present, in case labels differ from expectation

Uses `std_name` rather than `targname` — one physical star carries exactly one
`std_name`, whereas `BD+28 4211` alone appears as seven different `targname`
strings.

---



## 3. `download_std_red.py` — select and download



### Options

**Choosing frames**


| Option                    | What it does                                                                                          |
| ------------------------- | ----------------------------------------------------------------------------------------------------- |
| `--stars NAME ...`        | by `std_name`, e.g. `feige110 g191b2b`. Validated — an unknown name lists what is available           |
| `--gratings G ...`        | red gratings, e.g. `RL RM1 RH2`. Not validated: a typo silently yields "nothing left after filtering" |
| `--ifu IFU ...`           | `Large` / `Medium` / `Small`. Not validated either                                                    |
| `--nights YYYY-MM-DD ...` | Keck nights. Validated. **Not the UT date in the filename** — see Gotchas                             |
| `--public-only`           | drop frames still inside their proprietary period instead of logging in for them                      |


**Calibrations**


| Option               | What it does                                                                                                                                                                                                                               |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--same-night`       | the main mode. Keep only frames whose own night supplies their calibrations, and fetch just those. See below                                                                                                                               |
| `--require TYPE ...` | which lamp types must come from the frame's own night. Default: all five of `arclamp contbars flatlamp domeflat twiflat`. `bias` **is not a choice** — it is always required and comes from a second query per frame's night (SHARED bias) |
| `--keep-blue-calibs` | keep KB calibrations too. Default is red channel only                                                                                                                                                                                      |
| `--calib`            | the blunt alternative: hand `calibfile=1` to pykoa, which fetches **every date KOA associates**, unfiltered. Conflicts with `--same-night`                                                                                                 |


**Where things go**


| Option               | Default                     | What it does                                                                                      |
| -------------------- | --------------------------- | ------------------------------------------------------------------------------------------------- |
| `--outdir DIR`       | `pyKOA/fits`                | science lands in `<outdir>/lev0`, the per-night tree in `<outdir>/by_night`                       |
| `--calib-outdir DIR` | `<outdir>/calib_same_night` | calibrations. Kept apart from `<outdir>/calib`, which is where an unfiltered `--calib` run writes |
| `--cookie PATH`      | `outputKC/koa_cookie.txt`   | KOA login cookie, for proprietary frames                                                          |


**Running**


| Option      | Default | What it does                                                                       |
| ----------- | ------- | ---------------------------------------------------------------------------------- |
| `--dry-run` |         | report what would be selected, then stop. Nothing is downloaded                    |
| `--jobs N`  | 8       | concurrent caliblist queries. Each takes ~11 s, and there is one per science frame |




### What `--same-night` does

Two independent gates, in order.

**Lamps, from KOA's association.** Which calibration belongs to which science
frame is the archive's own bookkeeping — it turns on `stateid`, a configuration
hash not present in the metadata table — so it must be asked for one frame at a
time, `nph-getCaliblist?instrument=kcwi&koaid=…`. Each reply is a small JSON
table (~55 KB) naming files, not the files themselves. Within a list, rows from
the frame's own night survive; the frame is dropped if what remains does not
cover `--require`. Lists are cached in `outputKC/caliblists` and in whatever an
earlier `--calib` run left in `<outdir>/calib`, so the ~11 s query is paid once.

**Bias, from a query for the night.** *Not* from the association, because the
frames the association offers are the wrong ones. `stateid` encodes IFU,
gratings and wavelengths — `STATENAM` reads like `Small_BL4500_RL8150` — but not
binning and not readout mode. So it will return a 2×2 FPCam test bias for a 1×1
Small science frame and call it a match. Over 218 feige110 frames the
association supplied at least one same-night red bias for 32 of them, and three
or more for 9.

What the observatory actually takes is a shared afternoon bias set, good for
every program that night at that binning. Those are found by asking for the
night — `koaimtyp='bias'`, `camera='RED'` — and keeping the ones whose `binning`
and `ampmode` match the science frame. All 218 frames have seven or more, on all
41 nights.

`ccdspeed` is queried and reported but never matched on: biases are taken at
speed 0 almost without exception while science frames run 0, 1 and 2, so
requiring it to agree empties the pool for every frame that is not speed 0.

### Yields

For feige110, 218 frames:


| `--require`                 | frames kept | calibs | of which bias |
| --------------------------- | ----------- | ------ | ------------- |
| all five lamps (default)    | 114         | 1126   | 252           |
| `arclamp contbars flatlamp` | 209         | 1731   | 378           |


`twiflat` is what binds the default down to 114 — it is often associated across
nights. Arc lamps carry the wavelength solution and flexure moves it overnight,
so same-night arcs genuinely matter; twilight flats are stable for weeks. Bias
drops nobody, in either case.

### One directory per night

PypeIt groups frames by instrument configuration — for `keck_kcrm` that is
`dispname`, `decker`, `binning`, `cenwave` — and no part of that key is a date.
Calibration groups are then "a simple grouping of frames with the same
configuration". Point it at one flat directory and it will happily build a
single master arc out of four nights of arcs, undoing the gate above. Measured
here: 7 of 27 configurations span more than one night, covering 56 of 114
frames.

So the frames are also laid out one night per directory:

```
pypeit_setup -s keck_kcrm -r fits/by_night/2023-12-17
```

Nothing is copied. Every file already sits in `lev0` or
`calib_same_night/lev0`, and these are **hardlinks** — another name for the same
inode — so a night costs no bytes and deleting one destroys nothing. That is
what makes the tree disposable, which it needs to be: it is derived from the
filtering, and changing `--require` repartitions it. Each run wipes and rebuilds
rather than patching, clearing only files matching a KOAID pattern. Anything
else you leave in a night directory stays.

---



## Output files

Everything except the FITS lands in `outputKC/`.


| File                            | Written by | Rows   | What it is                                                                                                                |
| ------------------------------- | ---------- | ------ | ------------------------------------------------------------------------------------------------------------------------- |
| `kcwi_red_cones_<starlist>.tbl` | fetch      | 1314   | raw archive reply, cached. Delete to re-query                                                                             |
| `std_red_matched.csv`           | fetch      | 1307   | **the table everything else reads.** 21 archive columns plus `std_name`, `sep_fov`, `in_fov`, `name_ok`                   |
| `std_red_manifest.tbl`          | fetch      | 1307   | the same frames as a KOA manifest, local labels stripped                                                                  |
| `download_manifest.tbl`         | download   | varies | the science frames this run will fetch                                                                                    |
| `night_bias.tbl`                | download   | 460    | raw bias query — every red bias on the selected nights, before matching. Audit what was rejected, not just what was taken |
| `calib_manifest.tbl`            | download   | varies | the calibrations this run will fetch, deduped. Only exists after a real download                                          |
| `same_night_pairs.csv`          | download   | 3715   | one row per science↔calibration link                                                                                      |
| `caliblists/*.json`             | download   | —      | cached association replies, one per science frame                                                                         |




### `same_night_pairs.csv`

The download list is deduped — a calibration shared by twelve frames is fetched
once — so **this is the only place the many-to-many pairing survives**.

```
sci_koaid,std_name,night,rgratnam,ifunam,calib_koaid,koaimtyp,calib_date_obs,calib_ut,calib_elaptime,source
```

`source` says which mechanism found the row, `caliblist` or `night`. The two are
not interchangeable — one is KOA's judgement about this frame, the other is ours
about the night — and once deduped into a single manifest nothing else on disk
records the difference. Every `bias` row is `night`-sourced.

### The FITS


| Directory                     | What                                                                    |
| ----------------------------- | ----------------------------------------------------------------------- |
| `fits/lev0/`                  | science frames                                                          |
| `fits/calib_same_night/lev0/` | their calibrations, lamps and biases together                           |
| `fits/by_night/<night>/`      | hardlinks to both, one directory per night — **point PypeIt here**      |
| `fits/calib/`                 | only from an unfiltered `--calib` run. Also holds cached caliblist JSON |


---



## Typical workflows

**See what exists, download nothing**

```bash
python fetch_std_red.py                 # once; ~9 s, writes std_red_matched.csv
cd pyKOA && python analyze_gratings.py  # which star covers which grating
```

**Check what a selection would give you**

```bash
python download_std_red.py --stars feige110 --same-night --dry-run
python download_std_red.py --stars feige110 --same-night \
       --require arclamp contbars flatlamp --dry-run
```

The dry run prints the gate outcome — kept, dropped and why, per type — then the
selection by star and grating. That is what it is for: seeing how many frames
survive the filtering.

**One night, ready to reduce**

```bash
python download_std_red.py --stars feige110 --same-night --nights 2023-12-17
pypeit_setup -s keck_kcrm -r fits/by_night/2023-12-17
```

**Everything for one standard**

```bash
python download_std_red.py --stars feige110 --same-night
```

Interrupt it freely. pykoa skips any koaid already on disk, so re-running
resumes, and the per-night tree is rebuilt from scratch each time.

---



## Gotchas

**A Keck night is not a UT date.** Keck is UT−10, so one night straddles two UT
dates. `keck_night()` shifts timestamps back 22 h — 10 to reach HST, 12 more to
move the boundary to local noon — so a night's afternoon calibrations, the
observing, and the morning wrap-up all carry one label. This is why
`--nights 2023-12-17` returns frames named `KR.20231218.*`.

KOA's own `diff_date` column disagrees with this on 354 of 22803 associations —
every one of them a calibration taken between 12:00 and 14:00 HST. `diff_date`
subtracts UT dates, so on all 354 it is the one that is wrong: it calls a bias
15 h away from its science frame "same date".

**Downloads never overwrite.** pykoa skips any file already present and deletes
nothing, so `lev0` and `calib_same_night/lev0` only grow. Change `--require` and
the frames the old rule chose stay behind while the manifests are rewritten to
describe the new one — so the script reports files present that this run's
manifest does not name. It reports rather than deletes: they may belong to
another selection.

`--gratings` **and** `--ifu` **are not validated.** A typo produces "nothing left
after filtering", not an error naming the bad value. `--stars` and `--nights`
do check.

**Frames flagged** `bad` **by the archive still pass the lamp gate.** They are never
*required*, but nothing excludes them either, so they can reach the download.
None appear in the feige110 selection.

**The drop counts in the gate report do not sum to the total dropped.** A frame
missing three types is counted under all three. They are per-type frame counts,
not a partition.