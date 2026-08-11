"""Find and download all KCWI *red-channel* frames for the selected standard stars.

Strategy: OR the starlist's cone conditions into a single ADQL WHERE clause so the
server does the spatial cut.  One request, ~9 s, ~600 KB back.  Local matching is
then only needed to label each frame with which standard it landed on.
"""

import os
import re

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from pykoa.koa import Koa

STARLIST = "kcwi_stds_selected.txt"
OUT = "pyKOA/outputKC"
FITSDIR = "pyKOA/fits"
# Server-side cone.  Must stay well wider than the FOV: it has to catch frames the
# observer took at a star's proper-motion-corrected position (up to ~30" away from
# the epoch-2000 coordinates in the starlist) before any local cut can see them.
MATCH_RADIUS = 0.05 * u.deg

# Local test: was the star actually inside the IFU?  Largest slicer is 33" x 20.4",
# so its half-diagonal is the furthest a pointing can be and still hold the star.
FOV_RADIUS = 19.4

os.makedirs(OUT, exist_ok=True)

# NB: do NOT filter on the `object` string.  Both OBJECT and TARGNAME are sticky --
# whichever the observer last set persists into following exposures -- so a frame
# labelled "DOME/Dark/240.0 sec." or "NGC2775sky" is often a real on-target
# standard exposure carrying a stale label.  Measured: 29 of 40 DOME/Dark frames
# here have an elaptime that contradicts the exposure time in their own label, and
# some (hd93521) show 80x the flux of a typical frame.  koaimtyp/imtype/caltype and
# the lamp state are identical for real and mislabelled frames, so they cannot
# separate them either.  RA/DEC tracks TARGNAME, so the cone cut plus raoff/decoff
# is the trustworthy selector.

# --- 1. the standards -------------------------------------------------------
PAT = re.compile(
    r'^(?P<name>\S+)\s+'
    r'(?P<ra>\d+\s+\d+\s+[\d.]+)\s+'
    r'(?P<dec>[+-]?\d+\s+\d+\s+[\d.]+)\s+'
    r'(?P<epoch>[\d.]+)'
    r'.*?pmra=\s*(?P<pmra>[-+\d.eE]+)'
    r'.*?pmdec=\s*(?P<pmdec>[-+\d.eE]+)'
)

ALIAS_IN_COMMENT = re.compile(r'=\s*([A-Za-z0-9 +-]+?)\s*$')

# Abbreviations no rule recovers; extend as you meet them.
EXTRA_ALIAS = {"hil600": ["hilt 600"]}          # Hiltner 600


def norm(s):
    """Reduce a star name to a comparable form.

    The starlist spells the sign of the declination zone as a letter -- 'd' for
    '+', 'm' for '-' -- so bd02d3375 and an observer's BD023375 or BD+02 3375 all
    collapse to bd023375.
    """
    s = re.sub(r'[^a-z0-9]', '', str(s).lower())
    return re.sub(r'(?<=\d)[dm](?=\d)', '', s)


rows = []
aliases = {}
with open(STARLIST) as f:
    for i, line in enumerate(f, 1):
        body, _, comment = line.partition('#')
        body = body.strip()
        if not body:
            continue
        m = PAT.match(body)
        if not m:
            raise ValueError(f"line {i} unparsed: {body!r}")
        d = m.groupdict()
        rows.append(dict(name=d['name'], ra=d['ra'], dec=d['dec'],
                         pmra=float(d['pmra']), pmdec=float(d['pmdec'])))
        # the starlist records some cross-identifications in the comment: "... =GD 279"
        known = {norm(d['name'])} | {norm(a) for a in EXTRA_ALIAS.get(d['name'], [])}
        a = ALIAS_IN_COMMENT.search(comment.strip())
        if a:
            known.add(norm(a.group(1)))
        aliases[d['name']] = known

stds = pd.DataFrame(rows)
std_c = SkyCoord(stds.ra.values, stds.dec.values, unit=(u.hourangle, u.deg))

# --- 2. one server-side query: every red science frame inside any std's cone --
# koaid, instrume and filehand are mandatory if you want Koa.download() to work.
# elaptime is the exposure time -- `exptime` is entirely null for KCWI.
# raoff/decoff show whether the telescope was offset from the target.
# The rest is instrument configuration, needed downstream to pair a standard with
# science frames taken in the same setup.  koaimtyp and camera are not selected:
# the WHERE clause already pins them to 'object' and 'RED'.
COLS = """koaid, instrume, filehand, targname,
          rgratnam, rcwave, rgrangle, rnasnam, ifunam, binning, ampmode,
          elaptime, ra, dec, raoff, decoff, utdatetime,
          progid, progpi, propint, filesize_mb"""

cones = " or ".join(
    f"contains(point('J2000',ra,dec),"
    f"circle('J2000',{r:.6f},{d:.6f},{MATCH_RADIUS.to_value(u.deg)}))=1"
    for r, d in zip(std_c.ra.deg, std_c.dec.deg))

# Name the cache after the starlist: the nearest-neighbour labelling below is only
# valid while the cones that produced this table match the stars we label against.
index_tbl = f"{OUT}/kcwi_red_cones_{os.path.splitext(os.path.basename(STARLIST))[0]}.tbl"
if not os.path.exists(index_tbl):
    Koa.query_adql(
        f"select {COLS} from koa_kcwi "
        f"where koaimtyp = 'object' and camera = 'RED' and ({cones}) "
        f"order by utdatetime",
        index_tbl, overwrite=True, format="ipac")

red = Table.read(index_tbl, format="ascii.ipac").to_pandas()
# IPAC is fixed-width, so names and values arrive padded to the column width.
# Not select_dtypes("object"): pandas 3 gives strings their own 'str' dtype and
# will stop returning them for "object", leaving the padding in place.
red.columns = red.columns.str.strip()
for c in red.columns:
    if red[c].dtype == object or red[c].dtype.name == "str":
        red[c] = red[c].str.strip()

red_c = SkyCoord(red.ra.values * u.deg, red.dec.values * u.deg)

# --- 3. label each frame with which standard it hit ------------------------
idx, sep, _ = red_c.match_to_catalog_sky(std_c)
matched = red.copy()
matched["std_name"] = stds.name.values[idx]

# --- was the star inside the IFU? -------------------------------------------
# Observers were inconsistent about correcting for proper motion at the telescope:
# some pointed at the epoch-2000 coordinates, others at the star's true position
# that night.  Neither reference alone is right, so accept whichever is closer --
# testing only against the stale coordinates would reject the frames that were
# acquired correctly (wd1327m083 sits 30" from its epoch-2000 position).
obs_year = pd.to_datetime(matched.utdatetime, format="mixed")
dt = (obs_year.dt.year + obs_year.dt.dayofyear / 365.25).values - 2000.0
p = stds.set_index("name").loc[matched.std_name]
now_c = SkyCoord(
    (std_c.ra.deg[idx] + p.pmra.values * 15 * dt / 3600) * u.deg,   # pmra: sec of time/yr
    (std_c.dec.deg[idx] + p.pmdec.values * dt / 3600) * u.deg)
matched["sep_fov"] = np.minimum(sep.arcsec, red_c.separation(now_c).arcsec)
matched["in_fov"] = matched.sep_fov <= FOV_RADIUS

# --- does the observer's TARGNAME name this standard? ------------------------
matched["name_ok"] = [
    any(a in norm(t) or norm(t) in a for a in aliases[s])
    for s, t in zip(matched.std_name, matched.targname)]

print(f"{len(matched)} red frames on {matched.std_name.nunique()} standards")
print(matched.groupby("std_name").size().sort_values(ascending=False).head(20))
print(matched.groupby("rgratnam").size())
print(f"\nstar inside the {FOV_RADIUS}\" IFU field (in_fov): "
      f"{matched.in_fov.sum()} of {len(matched)}, "
      f"{matched[matched.in_fov].std_name.nunique()} standards")
print(f"TARGNAME names the matched standard (name_ok): "
      f"{matched.name_ok.sum()} of {len(matched)}")
print(f"telescope offset from target (|raoff|+|decoff| > 0): "
      f"{((matched.raoff.abs() + matched.decoff.abs()) > 0).sum()} frames")
if (~matched.name_ok).any():
    print("\nTARGNAME not recognised (add to EXTRA_ALIAS if these are the same star):")
    print(matched[~matched.name_ok].groupby(["std_name", "targname"])
          .size().rename("frames").to_string())

# --- keep only frames that really are the standard ---------------------------
# Precision over recall: there are plenty of frames per star, so a doubtful one is
# cheaper to drop than to carry.  A commanded offset needs no separate rule --
# ra/dec already includes it, so in_fov catches a nod that took the star off-slice.
keep = matched.name_ok & matched.in_fov
dropped = matched[~keep]
matched = matched[keep]

print(f"\ndropped {len(dropped)} frames:")
for lab, cond in [("TARGNAME mismatch", ~dropped.name_ok),
                  (f"star outside the {FOV_RADIUS}\" field", ~dropped.in_fov)]:
    print(f"  {cond.sum():>4}  {lab}")
print(f"kept {len(matched)} frames on {matched.std_name.nunique()} standards, "
      f"{matched.filesize_mb.sum() / 1024:.1f} GB")

thin = matched.groupby("std_name").size().sort_values()
print(f"\nthinnest coverage (check these before relying on them):")
print(thin.head(5).to_string())

# --- 4. write a download manifest in the format Koa.download() expects ------
# The manifest carries the KOA columns only -- the labels added above are local and
# have no meaning to the archive.  They stay in the CSV, which is the file to read
# when you want to know which standard a koaid belongs to.
LABELS = ["std_name", "sep_fov", "in_fov", "name_ok"]
manifest = f"{OUT}/std_red_manifest.tbl"
Table.from_pandas(matched.drop(columns=LABELS)).write(
    manifest, format="ascii.ipac", overwrite=True)
matched.to_csv(f"{OUT}/std_red_matched.csv", index=False)

# --- 5. download ------------------------------------------------------------
# propint = proprietary period in months; rows still inside it need Koa.login().
# calibfile=1 also drags down the associated biases/flats/arcs for each frame.
# Koa.download(manifest, "ipac", FITSDIR, calibfile=0)
