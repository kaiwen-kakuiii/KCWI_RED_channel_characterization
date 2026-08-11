"""Download the KCWI red-channel standard-star frames selected by fetch_std_red.py.

Reads std_red_matched.csv -- the labelled table, the one that knows which standard
each frame belongs to -- lets you cut it down to the stars/gratings/IFU you actually
want, writes that subset as a KOA manifest, and hands it to Koa.download().

    python download_std_red.py --dry-run                      # 1307 frames, 18.8 GB
    python download_std_red.py --stars feige110 g191b2b       # by standard
    python download_std_red.py --gratings RL RM1 --ifu Large  # by setup
    python download_std_red.py --public-only                  # skip the 3 proprietary
    python download_std_red.py --calib                        # + biases/flats/arcs

Re-running is safe and resumes: pykoa skips any koaid already on disk, so an
interrupted download continues where it stopped.

--same-night
------------
Keeps only frames whose calibrations come from their own observing night, and
fetches just those calibrations:

    python download_std_red.py --stars feige110 --same-night --dry-run
    python download_std_red.py --stars feige110 --same-night --require arclamp contbars flatlamp

Which calibration belongs to which science frame is KOA's own bookkeeping, not
something derivable from the columns queried here -- it turns on `stateid`, a hash
of the instrument configuration that is not even in the metadata table.  So the
association is asked for one koaid at a time, from

    nph-getCaliblist?instrument=kcwi&koaid=<koaid>

which returns a small JSON table (~55 KB) naming the associated files.  Those lists
are all that gets fetched up front; the FITS are downloaded only after filtering.
For feige110 that is 12 MB of JSON to decide 90 downloads instead of 3033.  Lists
are cached, in outputKC/caliblists and in whatever an earlier --calib run left in
<outdir>/calib, so the queries are paid once (~11 s each, run --jobs at a time).

Two filters come out of one pass over those lists:

  * within a list, keep the rows from the science frame's own night, red channel
    only unless --keep-blue-calibs;
  * drop the science frame entirely if what is left does not cover --require.

"Same night" is a Keck night, not a UT date.  Keck is UT-10, so afternoon calibs
taken before 14:00 HST carry the previous UT date while belonging to the coming
night; keck_night() shifts by 22 h to put the boundary at local noon.  Measured
over 22803 associations, KOA's own diff_date column disagrees with this on 354 --
every one of them taken between 12:00 and 14:00 HST, the boundary this exists to
place.  diff_date subtracts UT dates, so on all 354 it is the one that is wrong.

Requiring all five lamp types same-night is strict, because twiflat is often
associated across nights -- feige110 keeps 114 of 218 frames under the full set
and 209 under arclamp+contbars+flatlamp.  Both are defensible: arc lamps carry the
wavelength solution and flexure moves it overnight, while twilight flats are
stable for weeks.  --require is how you choose.

bias
----
Bias is not in --require and does not come from the caliblist, because the frames
the caliblist offers are the wrong ones.  KOA associates by stateid, which encodes
the IFU, the gratings and their wavelengths -- STATENAM reads like
'Small_BL4500_RL8150' -- but not the binning and not the readout mode.  So the
association returns a 2x2 FPCam test bias for a 1x1 Small science frame and calls
it a match; for the 2023-12-17 frames here that lone test bias was the entire
same-night bias set.  Over 218 feige110 frames the association supplied at least
one same-night red bias for 32 of them, and three or more for 9.

What the observatory actually takes is a shared afternoon bias set, good for every
program that night at that binning.  Those are found by asking for the night
rather than for the frame -- koaimtyp='bias', camera='RED' -- and keeping the ones
whose binning and ampmode match the science frame.  All 218 frames have seven or
more, on all 41 nights.  Bias therefore never drops a frame here, but the check
stays: a night with no matching bias is a night that cannot be reduced.

ccdspeed is queried and reported but never matched on.  Biases are taken at speed
0 almost without exception while science frames run 0, 1 and 2, so requiring it to
agree empties the pool for every frame that is not speed 0 -- 34 of 485 on three
sampled nights, each of which has seven perfectly good biases.

one directory per night
-----------------------
PypeIt groups frames by instrument configuration -- for keck_kcrm dispname,
decker, binning and cenwave -- and no part of that key is a date.  Calibration
groups are then "a simple grouping of frames with the same configuration", so a
single flat directory invites it to build one master arc out of four nights of
arcs, undoing the gate above.  Here 7 of 27 configurations span more than one
night, covering 56 of 114 frames.

So the frames are also laid out one night per directory, under <outdir>/by_night,
and pypeit_setup is pointed at a single night:

    pypeit_setup -s keck_kcrm -r fits/by_night/2023-12-17

Nothing is copied.  Every file already sits in lev0 or calib_same_night/lev0, and
these are hardlinks -- another name for the same inode -- so a night costs no
bytes and removing one destroys nothing.  That is what makes the tree disposable,
which it needs to be: it is derived from the filtering, and changing --require or
the bias match repartitions it, so each run wipes and rebuilds rather than
patching.  Only files this script could have created are cleared; anything else
left in a night directory stays.

    python download_std_red.py --stars feige110 --same-night --nights 2023-12-17

--nights takes the Keck night, not the UT date in the koaid.  The frames above are
KR.20231218.*, and they belong to the night of 2023-12-17.

Outputs, beyond the FITS themselves:

    outputKC/calib_manifest.tbl     the same-night calibrations, deduped
    outputKC/night_bias.tbl         the raw bias query, one row per bias on file
    outputKC/same_night_pairs.csv   one row per science/calibration link

The download list is deduped -- a calibration shared by several frames is fetched
once -- so the pairs CSV is the only place the many-to-many pairing survives.  Its
`source` column says which mechanism found each row, caliblist or night query.
"""

import argparse
import collections
import json
import os
import re
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

import pandas as pd
from astropy.table import Table
from dateutil.relativedelta import relativedelta
from pykoa.koa import Koa

# Paths resolve against this file, not the shell's cwd -- fetch_std_red.py is run
# from the repo root and analyze_gratings.py from pyKOA/, and this script should
# work from either.
HERE = os.path.dirname(os.path.abspath(__file__))
MATCHED_CSV = os.path.join(HERE, "outputKC", "std_red_matched.csv")
FITSDIR = os.path.join(HERE, "fits")
COOKIE = os.path.join(HERE, "outputKC", "koa_cookie.txt")

# Columns added locally by fetch_std_red.py; the archive has no use for them.
LABELS = ["std_name", "sep_fov", "in_fov", "name_ok"]

# KOA decides which calibrations belong to a science frame; the association is not
# derivable from the metadata we hold, so it has to be asked for one koaid at a
# time.  The reply is a small JSON table (~55 KB), not the FITS files themselves.
CALIBLIST_URL = ("https://koa.ipac.caltech.edu/cgi-bin/KoaAPI/nph-getCaliblist?"
                 "instrument=kcwi&koaid=")
CALIB_CACHE = os.path.join(HERE, "outputKC", "caliblists")
# Lamp calibrations only.  Bias is absent deliberately -- see the module docstring:
# the association returns test biases, so they are asked for by night instead.
LAMP_TYPES = ["arclamp", "contbars", "flatlamp", "domeflat", "twiflat"]

# binning and ampmode are the match keys; ccdspeed is carried for the report only.
BIAS_COLS = ("koaid, instrume, filehand, koaimtyp, camera, binning, ampmode, "
             "ccdspeed, ifunam, elaptime, utdatetime, progid, propint, filesize_mb")

# What the per-night tree is allowed to delete when it rebuilds: files it could
# have created itself, and nothing else.
KOAID_RE = re.compile(r"^K[RB]\.\d{8}\.\d+\.\d+\.fits$")


def keck_night(dt):
    """Label the Keck observing night a UT timestamp belongs to.

    Not the UT date: Keck is UT-10, so afternoon calibrations taken before 14:00
    HST fall on the previous UT date while belonging to the coming night.  Shift
    by 22 h -- 10 to reach HST, 12 more to move the boundary to local noon -- and
    one night gets one label.  Measured on 7922 associations, KOA's own diff_date
    field disagrees with this on 18 of them.
    """
    return (dt - timedelta(hours=22)).date()


def caliblist(koaid, cache_dirs):
    """KOA's calibration association for one science frame, cached on disk.

    Cached under either name pykoa or this script uses, so an earlier --calib run
    is reused rather than re-fetched.  Returns None when the archive cannot be
    reached, which is kept distinct from "no calibrations": one is a network
    failure to retry, the other is a real property of the frame.

    An unreadable cache file counts as a miss, not an error.  Interrupting a run
    mid-write used to leave a truncated JSON behind, and since the read was
    unguarded that file raised inside a worker thread and killed every later run
    too -- the one failure a re-run could not clear.  A bad copy is skipped rather
    than fetched around, so the other cache dir still answers when it can; only
    when none of them do does this go back to the archive.  Either way the write
    below repairs CALIB_CACHE, which is searched first.
    """
    stem = koaid[:-5] if koaid.endswith(".fits") else koaid
    for d in cache_dirs:
        p = os.path.join(d, f"{stem}.caliblist.json")
        if os.path.exists(p):
            try:
                with open(p) as f:
                    return json.load(f).get("table", [])
            except (OSError, ValueError):     # JSONDecodeError is a ValueError
                continue
    try:
        with urllib.request.urlopen(CALIBLIST_URL + koaid, timeout=120) as r:
            payload = json.load(r)
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
        return None
    # Write to a private name and rename: os.replace is atomic, so a reader either
    # sees the old file or the whole new one, never a half-written one.  The pid
    # keeps two concurrent runs of this script off each other's temporaries.
    os.makedirs(CALIB_CACHE, exist_ok=True)
    final = os.path.join(CALIB_CACHE, f"{stem}.caliblist.json")
    tmp = f"{final}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w") as f:
            json.dump(payload, f)
        os.replace(tmp, final)
    except OSError:
        try:
            os.unlink(tmp)
        except OSError:
            pass
    return payload.get("table", [])


def same_night(rows, night, keep_blue):
    """The calibrations from this frame's own night, red channel unless asked."""
    out = []
    for r in rows:
        if not keep_blue and not r["koaid"].startswith("KR"):
            continue
        t = datetime.strptime(f"{r['date_obs']} {r['ut'][:8]}", "%Y-%m-%d %H:%M:%S")
        if keck_night(t) == night:
            out.append(r)
    return out


def night_biases(nights, out_tbl):
    """Every red bias KOA holds on these Keck nights, whoever it was taken for.

    One query for all of them, not one per night: the reply is a few hundred rows
    and the round trip dominates.  A Keck night straddles two UT dates and koaid
    carries the UT date, so both are asked for and keck_night makes the real cut.
    Selecting on utdatetime would be the obvious move and does not work -- the
    archive is Oracle behind ADQL and rejects the string as a date (ORA-01843) --
    so the date lives in a `koaid like` pattern instead.
    """
    uts = sorted({(n + timedelta(days=k)).strftime("%Y%m%d")
                  for n in nights for k in (0, 1)})
    like = " or ".join(f"koaid like 'KR.{u}.%'" for u in uts)
    Koa.query_adql(f"select {BIAS_COLS} from koa_kcwi "
                   f"where koaimtyp = 'bias' and camera = 'RED' and ({like}) "
                   f"order by utdatetime",
                   out_tbl, overwrite=True, format="ipac")
    b = Table.read(out_tbl, format="ascii.ipac").to_pandas()
    b.columns = b.columns.str.strip()
    for c in b.columns:
        if b[c].dtype == object or b[c].dtype.name == "str":
            b[c] = b[c].str.strip()
    b["night"] = pd.to_datetime(b.utdatetime, format="mixed").apply(keck_night)
    return b


def bias_row(r):
    """A bias query row in the shape the caliblist uses, so the two can be mixed.

    The sources describe the same thing in different words: the query returns
    utdatetime where a caliblist returns date_obs and ut, and knows nothing of
    stateid.  Reconciling it once here is what lets the manifest, the dedup and
    the download stay unaware of where a row came from.
    """
    d, _, u = str(r.utdatetime).partition(" ")
    return {"koaid": r.koaid, "instrume": r.instrume, "filehand": r.filehand,
            "koaimtyp": "bias", "object": "", "targname": "", "date_obs": d,
            "ut": u or "00:00:00.00", "elaptime": r.elaptime, "camera": "RED",
            "rgratnam": "", "ifunam": r.ifunam, "stateid": "", "progid": r.progid,
            "semid": "", "propint": r.propint, "_source": "night"}


def load_matched(path):
    """Read the labelled CSV, undoing the fixed-width padding it inherits from IPAC.

    skipinitialspace is not optional: `binning` is written as the quoted field
    "2,2", and without it that comma splits the row and shifts every later column.
    """
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    # Not select_dtypes("object"): pandas 3 gives strings their own 'str' dtype and
    # will stop returning them for "object", at which point nothing would be
    # stripped and every name comparison below would silently miss.
    for c in df.columns:
        if df[c].dtype == object or df[c].dtype.name == "str":
            df[c] = df[c].str.strip()
    return df


def released(df):
    """True where the proprietary period has expired.

    propint is that period in months, counted from the observation.  Rows still
    inside it need a KOA account that owns the program; without a cookie the
    archive returns an error page instead of the FITS file.
    """
    obs = pd.to_datetime(df.utdatetime, format="mixed")
    now = datetime.now()
    return [o + relativedelta(months=int(m)) <= now
            for o, m in zip(obs, df.propint)]


ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--stars", nargs="+", metavar="NAME",
                help="standards to fetch, by std_name (default: all)")
ap.add_argument("--gratings", nargs="+", metavar="G",
                help="red gratings to keep, e.g. RL RM1 RH2 (default: all)")
ap.add_argument("--ifu", nargs="+", metavar="IFU", help="Large / Medium / Small")
ap.add_argument("--nights", nargs="+", metavar="YYYY-MM-DD",
                help="Keck nights to keep -- the night the frame belongs to, not "
                     "its UT date (they differ after 14:00 HST)")
ap.add_argument("--outdir", default=FITSDIR, help=f"default: {FITSDIR}")
ap.add_argument("--calib", action="store_true",
                help="also fetch the associated biases/flats/arcs (much larger) -- "
                     "every date KOA associates, not just the frame's own night")
ap.add_argument("--same-night", action="store_true",
                help="keep only frames whose own night supplies the calibrations "
                     "listed in --require, and fetch just those calibrations")
ap.add_argument("--require", nargs="+", metavar="TYPE", default=LAMP_TYPES,
                choices=LAMP_TYPES,
                help="lamp calib types that must come from the science frame's "
                     f"own night (default: all of {' '.join(LAMP_TYPES)}).  bias "
                     "is not listed: it is always required, and comes from a "
                     "query for the night rather than from the association")
ap.add_argument("--keep-blue-calibs", action="store_true",
                help="keep KB calibrations too (default: red channel only)")
ap.add_argument("--jobs", type=int, default=8, metavar="N",
                help="concurrent caliblist queries (default 8; ~11 s each)")
ap.add_argument("--calib-outdir", metavar="DIR",
                help="default: <outdir>/calib_same_night")
ap.add_argument("--public-only", action="store_true",
                help="drop still-proprietary rows instead of logging in")
ap.add_argument("--cookie", default=COOKIE, help="KOA cookie file path")
ap.add_argument("--dry-run", action="store_true",
                help="report what would be downloaded, then stop")
args = ap.parse_args()

if not os.path.exists(MATCHED_CSV):
    sys.exit(f"{MATCHED_CSV} not found -- run fetch_std_red.py first")

df = load_matched(MATCHED_CSV)
n_all = len(df)

# --- 1. cut down to what was asked for ---------------------------------------
if args.stars:
    want = {s.lower() for s in args.stars}
    unknown = want - set(df.std_name.str.lower())
    if unknown:
        sys.exit(f"no frames for: {', '.join(sorted(unknown))}\n"
                 f"available: {', '.join(sorted(df.std_name.unique()))}")
    df = df[df.std_name.str.lower().isin(want)]
if args.gratings:
    df = df[df.rgratnam.str.upper().isin({g.upper() for g in args.gratings})]
if args.ifu:
    df = df[df.ifunam.str.lower().isin({i.lower() for i in args.ifu})]
if args.nights:
    # Matched against the Keck night, not the UT date in the koaid: a frame taken
    # at 02:00 UT on the 18th belongs to the night of the 17th, and asking for
    # "2023-12-18" because that is what the filename says would return nothing.
    try:
        want = {datetime.strptime(n, "%Y-%m-%d").date() for n in args.nights}
    except ValueError as e:
        sys.exit(f"--nights wants YYYY-MM-DD: {e}")
    got = pd.to_datetime(df.utdatetime, format="mixed").apply(keck_night)
    unknown = want - set(got)
    if unknown:
        avail = sorted(set(got))
        sys.exit(f"no frames on: {', '.join(str(n) for n in sorted(unknown))}\n"
                 f"{len(avail)} nights available, "
                 f"{avail[0]} to {avail[-1]}")
    df = df[got.isin(want)]

if df.empty:
    sys.exit("nothing left after filtering")

# --- 2. proprietary rows ------------------------------------------------------
is_public = pd.Series(released(df), index=df.index)
n_prop = int((~is_public).sum())
n_dropped = 0
if args.public_only:
    df = df[is_public]
    n_prop, n_dropped = 0, n_prop      # they are gone, not pending a login
    if df.empty:
        sys.exit("every selected frame is still proprietary")

# --- 2b. same-night calibrations ----------------------------------------------
# Ask KOA which calibrations belong to each frame, keep the ones from that frame's
# own night, and drop the frame if they do not cover --require.  This runs before
# anything is sized or downloaded, so the gate decides what the rest of the script
# even considers.
pairs = {}                       # science koaid -> its same-night calib rows
if args.same_night:
    if args.calib:
        sys.exit("--same-night and --calib conflict: calibfile=1 fetches every "
                 "date KOA associates, which undoes the gate")
    required = set(args.require)
    nights = dict(zip(df.koaid,
                      pd.to_datetime(df.utdatetime, format="mixed").apply(keck_night)))
    cache_dirs = [CALIB_CACHE, os.path.join(args.outdir, "calib")]
    cached = sum(any(os.path.exists(os.path.join(d, f"{k[:-5]}.caliblist.json"))
                     for d in cache_dirs) for k in df.koaid)
    print(f"caliblists for {len(df)} frames: {cached} cached, "
          f"{len(df) - cached} to fetch (~{(len(df) - cached) * 11 / args.jobs / 60:.0f} min)")

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        lists = dict(zip(df.koaid,
                         ex.map(lambda k: caliblist(k, cache_dirs), df.koaid)))

    failed = [k for k, v in lists.items() if v is None]
    short = collections.Counter()
    keep = []
    for k, rows in lists.items():
        if rows is None:
            continue
        # Bias rows are dropped here, not merely ignored: they are the association's
        # test biases, and leaving them in the manifest would download the very
        # frames the night query exists to replace.  _source records where a row
        # came from, which koaimtyp cannot -- both mechanisms yield biases.
        mine = [dict(r, _source="caliblist")
                for r in same_night(rows, nights[k], args.keep_blue_calibs)
                if r["koaimtyp"] != "bias"]
        missing = required - {r["koaimtyp"] for r in mine}
        if missing:
            short.update(missing)
        else:
            keep.append(k)
            pairs[k] = mine

    print(f"\nsame-night gate (require: {' '.join(sorted(required))}"
          f"{'' if args.keep_blue_calibs else ', red only'})")
    print(f"   kept    {len(keep):>5} frames")
    print(f"   dropped {len(lists) - len(keep) - len(failed):>5} missing types "
          f"-- {dict(short.most_common())}")
    if failed:
        print(f"   dropped {len(failed):>5} caliblist fetch failed (re-run to retry)")
    if not keep:
        sys.exit("no frame has same-night calibrations covering --require")

    # --- bias, from the night rather than from the association ----------------
    # Runs after the lamp gate, on the frames that survived it: that is fewer
    # nights to ask about, and a bias could never rescue a frame the lamps have
    # already dropped.  One query covers every remaining night at once.
    meta_pre = df.set_index("koaid")
    want_nights = sorted({nights[k] for k in keep})
    bias_tbl = os.path.join(HERE, "outputKC", "night_bias.tbl")
    print(f"\nbias: querying KOA for every red bias on {len(want_nights)} nights ...")
    try:
        pool = night_biases(want_nights, bias_tbl)
    except Exception as e:
        # Unlike a caliblist, this is one request covering every frame.  There is
        # no partial result to carry on with, and continuing would silently mean
        # "no biases", so it stops here rather than quietly reducing the yield.
        sys.exit(f"bias query failed: {type(e).__name__}: {e}\n"
                 f"nothing was downloaded -- re-run to retry")

    nobias = collections.Counter()
    keep2 = []
    for k in keep:
        s = meta_pre.loc[k]
        m = pool[(pool.night == nights[k]) & (pool.binning == s.binning)
                 & (pool.ampmode == s.ampmode)]
        if m.empty:
            nobias[f"{s.binning} {s.ampmode}"] += 1
        else:
            keep2.append(k)
            pairs[k] = pairs[k] + [bias_row(r) for r in m.itertuples()]
    print(f"   {len(pool)} biases on those nights; "
          f"matched on binning+ampmode for {len(keep2)} of {len(keep)} frames")
    if nobias:
        print(f"   dropped {sum(nobias.values())} frames whose night has no bias "
              f"at their setup -- {dict(nobias.most_common())}")
    if not keep2:
        sys.exit("no frame has a same-night bias matching its binning and ampmode")
    keep = keep2
    df = df[df.koaid.isin(keep)]
    # Recount: n_prop was taken before the gate, and the gate typically drops most
    # of the selection, so the stale count would report -- and log in for -- frames
    # that are no longer in the manifest.  is_public is indexed on the pre-gate
    # frame and the gate only removes rows, so the lookup lines up.
    n_prop = int((~is_public[df.index]).sum())

# --- 3. report ----------------------------------------------------------------
# pykoa writes level-0 files to <outdir>/lev0 and skips any koaid already there.
lev0 = os.path.join(args.outdir, "lev0")
have = {f for f in os.listdir(lev0)} if os.path.isdir(lev0) else set()
todo = df[~df.koaid.isin(have)]

print(f"selected {len(df)} of {n_all} frames on {df.std_name.nunique()} standards")
print(df.groupby(["std_name", "rgratnam"]).size().rename("frames").to_string())
print(f"\nalready on disk: {len(df) - len(todo)}   to download: {len(todo)}"
      f"   ({todo.filesize_mb.sum() / 1024:.1f} GB)")
if n_prop:
    print(f"still proprietary: {n_prop} frames -- these need a KOA login "
          f"(or --public-only to skip them)")
if n_dropped:
    print(f"dropped {n_dropped} still-proprietary frames (--public-only)")
if args.calib:
    print("calibrations requested: expect several times the size above")

# One row per science-frame/calibration link.  The download list below is deduped
# -- a calib shared by several frames is fetched once -- so this file is the only
# place the many-to-many pairing survives.
calib_rows, pair_csv = {}, None
if args.same_night:
    pair_csv = os.path.join(HERE, "outputKC", "same_night_pairs.csv")
    meta = df.set_index("koaid")
    recs = []
    for sci, rows in pairs.items():
        if sci not in meta.index:
            continue
        for r in rows:
            calib_rows[r["koaid"]] = r
            # Which mechanism found this row.  The two are not interchangeable --
            # a caliblist row is KOA's judgement about this frame, a night row is
            # ours about the night -- and once deduped into one manifest nothing
            # else on disk records the difference.
            recs.append(dict(sci_koaid=sci, std_name=meta.std_name[sci],
                             night=nights[sci], rgratnam=meta.rgratnam[sci],
                             ifunam=meta.ifunam[sci], calib_koaid=r["koaid"],
                             koaimtyp=r["koaimtyp"], calib_date_obs=r["date_obs"],
                             calib_ut=r["ut"], calib_elaptime=r["elaptime"],
                             source=r["_source"]))
    pd.DataFrame(recs).to_csv(pair_csv, index=False)
    by_type = collections.Counter(r["koaimtyp"] for r in calib_rows.values())
    print(f"\nsame-night calibrations: {len(calib_rows)} distinct files "
          f"over {len(set(nights[k] for k in pairs))} nights, "
          f"{len(recs)} pair rows -> {os.path.relpath(pair_csv, HERE)}")
    print(f"   {dict(by_type.most_common())}")

if args.dry_run:
    print("\ndry run, stopping here.")
    sys.exit()
# Science frames all on disk does not mean there is nothing left to do: under
# --same-night the calibrations are a second download, and a re-run after the
# science finished is exactly when they are what remains.
if todo.empty and not calib_rows:
    print("\nnothing to do.")
    sys.exit()

# --- 4. login and the options shared by both downloads ------------------------
cookie = None
if n_prop and not args.public_only:
    cookie = args.cookie
    if not os.path.exists(cookie):
        print(f"\nKOA login for the {n_prop} proprietary frames "
              f"(cookie will be saved to {cookie}):")
        Koa.login(cookie)          # prompts for userid/password
        if not os.path.exists(cookie):
            sys.exit("login failed -- no cookie written")

kw = dict(lev0file=1, calibfile=1 if args.calib else 0)
if cookie:
    kw["cookiepath"] = cookie
os.makedirs(args.outdir, exist_ok=True)

# --- 5. the science frames ----------------------------------------------------
# Koa.download() needs instrume + filehand, and koaid too when it has to resolve
# each frame's calibrations.  The local label columns would only confuse it.
if todo.empty:
    print(f"\nall {len(df)} science frames already in {lev0}")
else:
    manifest = os.path.join(HERE, "outputKC", "download_manifest.tbl")
    Table.from_pandas(todo.drop(columns=[c for c in LABELS if c in todo])).write(
        manifest, format="ascii.ipac", overwrite=True)

    print(f"\ndownloading {len(todo)} frames into {lev0} ...")
    Koa.download(manifest, "ipac", args.outdir, **kw)

    now = {f for f in os.listdir(lev0)} if os.path.isdir(lev0) else set()
    got = todo.koaid.isin(now).sum()
    print(f"\n{got} of {len(todo)} files present in {lev0}")
    if got < len(todo):
        print("re-run to retry the rest -- finished files are skipped")

# --- 6. the same-night calibrations -------------------------------------------
# Lamps and biases go down together, in one manifest and one directory, even
# though they were found by different means -- the pairs CSV's `source` column is
# what records which.  Kept apart from <outdir>/calib, though: an earlier
# unfiltered --calib run writes every associated date there, and mixing the two
# would leave no way to tell a gated calibration from one that merely came close
# in time.
if calib_rows:
    cols = ["koaid", "instrume", "filehand", "koaimtyp", "object", "targname",
            "date_obs", "ut", "elaptime", "camera", "rgratnam", "ifunam",
            "stateid", "progid", "semid", "propint"]
    ctab = pd.DataFrame(list(calib_rows.values()))
    ctab = ctab[[c for c in cols if c in ctab]].astype(str)
    calib_manifest = os.path.join(HERE, "outputKC", "calib_manifest.tbl")
    Table.from_pandas(ctab).write(calib_manifest, format="ascii.ipac", overwrite=True)

    calib_out = args.calib_outdir or os.path.join(args.outdir, "calib_same_night")
    calib_lev0 = os.path.join(calib_out, "lev0")
    os.makedirs(calib_out, exist_ok=True)

    # Files here that this manifest does not name.  pykoa skips what already
    # exists and never deletes, so the directory only ever grows: change --require
    # or the bias match and the frames the old rule chose stay behind, while the
    # manifest and the pairs CSV are rewritten to describe the new one.  Nothing
    # downstream would notice -- pypeit_setup reads the directory, not the CSV --
    # so the drift is reported here rather than left to be discovered in a
    # reduction.  Reported, not deleted: they may be another selection's calibs.
    if os.path.isdir(calib_lev0):
        stale = sorted(set(os.listdir(calib_lev0)) - set(ctab.koaid))
        if stale:
            print(f"\n{len(stale)} file(s) in {calib_lev0} are not in this "
                  f"manifest, left by an earlier selection:")
            for s in stale[:5]:
                print(f"   {s}")
            if len(stale) > 5:
                print(f"   ... and {len(stale) - 5} more")
            print("   pypeit_setup reads the directory, so these would be picked "
                  "up too -- remove them if that is not what you want")

    print(f"\ndownloading {len(ctab)} same-night calibrations into {calib_lev0} ...")
    Koa.download(calib_manifest, "ipac", calib_out,
                 **{k: v for k, v in kw.items() if k != "calibfile"}, calibfile=0)

    on_disk = set(os.listdir(calib_lev0)) if os.path.isdir(calib_lev0) else set()
    print(f"{sum(k in on_disk for k in calib_rows)} of {len(calib_rows)} "
          f"calibrations present in {calib_lev0}")

    # --- 7. one directory per night -------------------------------------------
    # PypeIt groups frames by instrument configuration -- for keck_kcrm that is
    # dispname, decker, binning and cenwave -- and nothing in that key is a date.
    # Calibration groups are then "a simple grouping of frames with the same
    # configuration", so a flat directory lets it build one master arc from four
    # nights of arcs, which is exactly what the same-night gate exists to prevent.
    # Measured here: 7 of 27 configurations span more than one night, covering 56
    # of 114 frames.  Giving pypeit_setup one night at a time makes the merge
    # impossible rather than merely unlikely.
    #
    # Hardlinks, not copies: every file already sits in lev0 or calib_lev0, and a
    # hardlink is another name for the same inode, so a night costs no bytes and
    # unlinking one never destroys data.  That makes the tree disposable, which
    # matters because it is derived -- change --require or the bias match and the
    # nights repartition, so it is wiped and rebuilt rather than patched.
    tree = collections.defaultdict(set)
    for sci in df.koaid:
        n = str(nights[sci])
        tree[n].add(sci)
        tree[n].update(r["koaid"] for r in pairs.get(sci, []))

    by_night = os.path.join(args.outdir, "by_night")
    linked = absent = cleared = 0
    for n, koaids in sorted(tree.items()):
        d = os.path.join(by_night, n)
        os.makedirs(d, exist_ok=True)
        # Only KOAID-named regular files are cleared -- anything else you put in
        # here is yours and survives the rebuild.
        for f in os.listdir(d):
            p = os.path.join(d, f)
            if KOAID_RE.match(f) and os.path.isfile(p):
                os.unlink(p)
                cleared += 1
        for k in sorted(koaids):
            for src in (os.path.join(lev0, k), os.path.join(calib_lev0, k)):
                if os.path.exists(src):
                    try:
                        os.link(src, os.path.join(d, k))
                        linked += 1
                    except OSError as e:
                        # Cross-device is the one that actually happens, when
                        # --calib-outdir points at another filesystem.
                        print(f"   cannot link {k}: {e}")
                    break
            else:
                absent += 1

    print(f"\nper-night views in {by_night}: {linked} links over {len(tree)} nights"
          f"{f', {cleared} stale cleared' if cleared else ''}")
    if absent:
        print(f"   {absent} file(s) not linked -- not downloaded yet; "
              f"re-run to finish, the tree is rebuilt each time")
    print(f"   pypeit_setup -s keck_kcrm -r {os.path.join(by_night, sorted(tree)[0])}")
