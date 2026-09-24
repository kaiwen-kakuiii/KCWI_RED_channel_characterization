"""File logic shared by the throughput steps.

The scripts next to this file are the steps a person runs. This module is
the one implementation of the grating file, the .pypeit parameter block, the
exposure floor, and the visit split, so those cannot drift apart the way the
per-grating drivers did. It does not call PypeIt.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")

# A nod of 1.2 arcsec is one visit (RM2 2023-09-23 B). A repoint of 3.3
# arcsec is two visits (RH1 2023-10-15). 3 arcsec sits between those two.
POINTING_TOL_ARCSEC = 3.0
# Two sky windows that differ by at most this many percent of a slice share
# one reduction pass. Measured windows that really differ sit ~7 percent apart.
REGION_TOL_PCT = 3.0
# Instrument settings that are not in the grating file. Every grating uses them.
SPECTROGRAPH = "keck_kcrm"
SENS_ALGORITHM = "IR"
SENS_EXTR = "BOX"
TRIM_FRAC = 0.20
TRIM_PAD = 10
# A CDELT above this (degrees per spaxel) is the single-frame cube failure:
# the WCS is written at 1 degree and the extraction radius becomes huge.
CDELT_MAX_DEG = 0.01

REQUIRED = (
    "grating", "polyorder", "bridge", "lamps", "boxcar_arcsec",
    "exptime_floor_s", "min_frames", "max_airmass_spread", "count_rate_floor",
    "slice_span_A", "shift_tol_A", "nline_min", "templates",
)


class StepError(Exception):
    """A step stopped. The message is what the person running it should read."""


def pypeit_bin():
    return os.environ.get("PYPEIT_BIN", "/opt/miniconda3/envs/pypeit/bin")


def tool(name):
    return os.path.join(pypeit_bin(), name)


def log_path(grating, night):
    return os.path.join(night_dir(grating, night), "pipeline.log")


def next_step(script, grating, night, setup=None):
    """Print the next command using the Python that is running this step."""
    cmd = [sys.executable, os.path.join("pipeline", script),
           "--grating", grating, "--night", night]
    if setup is not None:
        cmd += ["--setup", setup]
    print("next:", " ".join(cmd), flush=True)


def run_logged(cmd, log, cwd=None, dry=False):
    shown = " ".join(cmd)
    print(f"$ {shown}", flush=True)
    if dry:
        return 0
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "a") as fh:
        fh.write(f"\n$ {shown}\n")
        fh.flush()
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    if proc.returncode:
        raise StepError(f"command failed (rc={proc.returncode}). See {log}")
    return proc.returncode


def read_pypeit(grating, night, setup, made_by):
    path = require_file(pypeit_path(grating, night, setup), made_by)
    with open(path) as fh:
        return path, fh.read()


def science_frames(text, cfg):
    """Active science rows after the exposure floor. Prints what it dropped."""
    rows = active(parse_rows(text), "science")
    if not rows:
        raise StepError("no science frames in this setup")
    gaps = calibration_gaps(parse_rows(text))
    if gaps:
        raise StepError(
            f"setup has science frames but no {', '.join(gaps)} frames. "
            f"Leave this setup; it is not a standard-star configuration.")
    kept, dropped = apply_exptime_floor(
        rows, cfg["exptime_floor_s"], cfg["min_frames"])
    if dropped:
        names = ", ".join(frame["filename"] for frame in dropped)
        print(f"exposure floor: {len(dropped)} frame(s) under "
              f"{float(cfg['exptime_floor_s']):.0f}s will not be reduced: {names}")
        print(f"{len(kept)} science frame(s) remain")
    elif any(frame["exptime"] < float(cfg["exptime_floor_s"]) for frame in rows):
        print(f"keeping frames under {float(cfg['exptime_floor_s']):.0f}s: "
              f"dropping them would leave fewer than {cfg['min_frames']}")
    return kept


def parse_yaml(text):
    """The subset the grating files use: maps, lists, scalars, comments."""
    root = {}
    # (indent, container). indent -1 is the root, so a top-level key stays.
    stack = [(-1, root)]
    pending = None  # (indent, parent_dict, key) whose value is the next block

    def flush_pending_as_dict():
        nonlocal pending
        if pending is None:
            return
        indent, parent, key = pending
        child = {}
        parent[key] = child
        stack.append((indent, child))
        pending = None

    for lineno, raw in enumerate(text.splitlines(), start=1):
        if "\t" in raw:
            raise StepError(f"line {lineno}: tabs are not allowed in a grating file")
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        line = raw.strip()

        if pending is not None and indent <= pending[0]:
            flush_pending_as_dict()
        while len(stack) > 1 and indent <= stack[-1][0]:
            stack.pop()

        if line.startswith("- "):
            if pending is not None and indent > pending[0]:
                indent_p, parent_p, key = pending
                child = []
                parent_p[key] = child
                stack.append((indent_p, child))
                pending = None
            parent = stack[-1][1]
            if not isinstance(parent, list):
                raise StepError(f"line {lineno}: list item is not under a list")
            parent.append(_scalar(line[2:].strip()))
            continue

        if pending is not None and indent > pending[0]:
            indent_p, parent_p, key_p = pending
            child = {}
            parent_p[key_p] = child
            stack.append((indent_p, child))
            pending = None

        parent = stack[-1][1]
        key, val = _split_kv(line, lineno)
        key = _unquote(key.strip())
        val = val.strip()
        if not isinstance(parent, dict):
            raise StepError(f"line {lineno}: expected a list item")
        if val == "":
            pending = (indent, parent, key)
        else:
            parent[key] = _scalar(val)
    flush_pending_as_dict()
    return root


def _split_kv(line, lineno):
    in_quote = None
    for i, ch in enumerate(line):
        if ch in "\"'" and in_quote is None:
            in_quote = ch
        elif ch == in_quote:
            in_quote = None
        elif ch == ":" and in_quote is None:
            return line[:i], line[i + 1:]
    raise StepError(f"line {lineno}: no key/value separator in {line!r}")


def _unquote(text):
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    return text


def _scalar(text):
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    if text in ("null", "~"):
        return None
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [_scalar(part) for part in inner.split(",")]
    try:
        if any(c in text for c in ".eE"):
            return float(text)
        return int(text)
    except ValueError:
        return text


def load_config(grating):
    path = os.path.join(CONFIG_DIR, f"{grating}.yaml")
    if not os.path.isfile(path):
        raise StepError(
            f"no grating file for {grating}. Expected {path}. "
            f"A new grating gets its own file; it does not borrow another grating's numbers.")
    with open(path) as fh:
        cfg = parse_yaml(fh.read())
    missing = [key for key in REQUIRED if key not in cfg]
    if missing:
        raise StepError(f"{path} is missing {', '.join(missing)}")
    if cfg["grating"] != grating:
        raise StepError(f"{path} says grating {cfg['grating']}")
    span = cfg["slice_span_A"]
    if not (isinstance(span, list) and len(span) == 2 and float(span[0]) < float(span[1])):
        raise StepError(f"{path}: slice_span_A must be [min, max] with min < max")
    if not isinstance(cfg["templates"], dict) or not cfg["templates"]:
        raise StepError(f"{path}: templates must list at least one central wavelength and slicer")
    for key, value in (cfg.get("exclude_regions") or {}).items():
        if not str(value).endswith(","):
            raise StepError(
                f"{path}: exclude_regions for {key} must end with a comma. "
                f"PypeIt reads the value as a list and walks a bare string one character at a time.")
    return cfg


def configured_lists():
    """Drop and flag tags declared in the grating files."""
    drops, flags = {}, {}
    if not os.path.isdir(CONFIG_DIR):
        return drops, flags
    for name in sorted(os.listdir(CONFIG_DIR)):
        if not name.endswith(".yaml"):
            continue
        cfg = load_config(name[:-5])
        for tag in cfg.get("drop") or []:
            drops.setdefault(tag, f"pipeline/config/{name}")
        for tag, reason in (cfg.get("flag") or {}).items():
            flags.setdefault(tag, reason)
    return drops, flags


def round_cenwave(value):
    """Nearest 10 Å. 8850, 8900 and 8950 must stay distinct."""
    return int(round(float(value) / 10.0) * 10)


def setting_key(cenwave, decker):
    return f"{int(cenwave)} {decker}"


def template_entry(cfg, cenwave, decker):
    key = setting_key(cenwave, decker)
    templates = cfg["templates"]
    if key not in templates:
        known = ", ".join(sorted(templates)) or "(none)"
        raise StepError(
            f"{key} is not in templates. Known settings: {known}. "
            f"Build a template for this central wavelength and slicer, add its "
            f"path to the grating file, then run this step again.")
    return templates[key]


def reid_value(grating, entry):
    """What `reid_arxiv` should say.

    A local file is an absolute path. PypeIt's own file is the bare filename:
    the cache copy has no extension, and writing that cache path makes PypeIt
    refuse to read it.
    """
    if entry == "shipped":
        return f"keck_kcrm_{grating}.fits"
    path = entry if os.path.isabs(entry) else os.path.join(ROOT, entry)
    if not os.path.isfile(path):
        raise StepError(f"template file does not exist: {path}")
    return os.path.abspath(path)


def lamps_value(name):
    """Absolute path PypeIt expects: the stem, without `_lines.dat`."""
    if name is None:
        return None
    stem = os.path.join(ROOT, "templates", str(name))
    listing = stem + "_lines.dat"
    if not os.path.isfile(listing):
        raise StepError(f"line list does not exist: {listing}")
    return os.path.abspath(stem)


def shift_tol(cfg, decker):
    raw = cfg["shift_tol_A"]
    if isinstance(raw, dict):
        if decker not in raw:
            raise StepError(
                f"no shift_tol_A for slicer {decker}. This grating's file gives "
                f"{', '.join(sorted(raw))}. Measure the blue-end spread on a "
                f"24/24 setup of this slicer and add it; do not reuse another slicer's number.")
        return float(raw[decker])
    return float(raw)


def setup_key(night, letter):
    return f"{night}/{letter}"


def optional_setup(cfg, name, night, letter):
    block = cfg.get(name) or {}
    return block.get(setup_key(night, letter))


def letter_of(setup):
    letter = setup.replace("keck_kcrm_", "")
    if not re.fullmatch(r"[A-Z]", letter):
        raise StepError(f"setup must be a letter such as B, not {setup!r}")
    return letter


def night_dir(grating, night):
    return os.path.join(ROOT, "reductions", grating, night)


def setup_dir(grating, night, setup):
    letter = letter_of(setup)
    return os.path.join(night_dir(grating, night), "pypeit_run", f"keck_kcrm_{letter}")


def pypeit_path(grating, night, setup):
    letter = letter_of(setup)
    return os.path.join(setup_dir(grating, night, setup), f"keck_kcrm_{letter}.pypeit")


def require_file(path, made_by):
    if not os.path.isfile(path):
        raise StepError(f"missing {path}. Run {made_by} first.")
    return path


def split_header(text):
    """Parameter block, then the setup and data blocks, which are kept verbatim."""
    match = re.search(r"\n# Setup\n", text)
    if match is None:
        match = re.search(r"\nsetup read\n", text)
    if match is None:
        raise StepError("no setup block in the .pypeit file")
    return text[:match.start()], text[match.start():]


def header_values(text):
    head, _ = split_header(text)
    found = {}
    for key in ("reid_arxiv", "lamps", "exclude_regions", "length_range",
                "user_regions", "method"):
        match = re.search(rf"^\s*{key}\s*=\s*(\S+)\s*$", head, re.M)
        found[key] = None if match is None else match.group(1)
    return found


def render_header(text, reid=None, lamps=None, exclude=None, length_range=None,
                  user_regions=None):
    """Replace the parameter block. The setup and data blocks are not rewritten."""
    head, body = split_header(text)
    banner = [line for line in head.splitlines()
              if line.startswith("# Auto-generated") or line.startswith("# UTC")]
    if not banner:
        banner = ["# Auto-generated PypeIt input file using PypeIt version: 2.0.1"]
    lines = banner + ["", "# User-defined execution parameters",
                      "[rdx]", "    spectrograph = keck_kcrm", ""]
    cal = []
    if exclude or length_range is not None:
        cal.append("    [[slitedges]]")
        if length_range is not None:
            cal.append(f"        length_range = {length_range}")
        if exclude:
            if not str(exclude).endswith(","):
                raise StepError("exclude_regions must end with a comma")
            cal.append(f"        exclude_regions = {exclude}")
    if reid or lamps:
        cal.append("    [[wavelengths]]")
        cal.append("        method = full_template")
        if reid:
            cal.append(f"        reid_arxiv = {reid}")
        if lamps:
            cal.append(f"        lamps = {lamps}")
    if cal:
        lines += ["[calibrations]"] + cal + [""]
    if user_regions:
        lines += ["[reduce]", "    [[skysub]]",
                  f"        user_regions = {user_regions}", ""]
    rendered = "\n".join(lines).rstrip("\n") + "\n"
    return rendered + body.lstrip("\n")


def setup_meta(text):
    match = re.search(r"^Setup\s+(\w+):\n(.*?)^setup end", text, re.M | re.S)
    if match is None:
        raise StepError("no Setup block")
    fields = {"letter": match.group(1)}
    for line in match.group(2).splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        fields[key.strip()] = value.strip()
    for key in ("cenwave", "decker", "dispname"):
        if key not in fields:
            raise StepError(f"Setup block has no {key}")
    return fields


def raw_dir(text):
    match = re.search(r"^\s*path\s+(\S+)\s*$", text, re.M)
    if match is None:
        raise StepError("data block has no path line")
    return match.group(1)


def parse_rows(text):
    """Rows of the data block. Commented rows are returned with commented=True."""
    match = re.search(r"^data read\n(.*)^data end", text, re.M | re.S)
    if match is None:
        raise StepError("no data block")
    columns = None
    rows = []
    for line in match.group(1).splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        commented = stripped.startswith("#")
        bare = stripped[1:].strip() if commented else stripped
        parts = [part.strip() for part in bare.split("|")]
        if columns is None:
            if parts and parts[0] == "filename":
                columns = parts
            continue
        if not parts or not parts[0].endswith(".fits"):
            continue
        row = dict(zip(columns, parts))
        row["commented"] = commented
        row["filename"] = parts[0]
        row["tokens"] = [tok for tok in row.get("frametype", "").split(",") if tok]
        try:
            row["exptime"] = float(row.get("exptime", "nan"))
        except ValueError:
            row["exptime"] = float("nan")
        try:
            row["ra"] = float(row.get("ra_off") or 0.0)
            row["dec"] = float(row.get("dec_off") or 0.0)
        except ValueError:
            row["ra"], row["dec"] = 0.0, 0.0
        row["star"] = _star(row.get("target", ""))
        rows.append(row)
    return rows


def _star(name):
    return re.sub(r"[^a-z0-9]+", "", (name or "").strip().lower()) or "std"


def star_name(name):
    """Target name used in cube tags: lowercase, letters and digits only."""
    return _star(name)


def active(rows, token):
    return [row for row in rows
            if not row["commented"] and token in row["tokens"]]


def calibration_gaps(rows):
    tokens = {tok for row in rows if not row["commented"] for tok in row["tokens"]}
    return [name for name in ("arc", "trace", "pixelflat", "align") if name not in tokens]


def apply_exptime_floor(frames, floor, min_frames):
    """Drop frames shorter than `floor` unless that would leave fewer than `min_frames`.

    Returns (kept, dropped). An empty dropped list means the floor yielded.
    """
    short = [frame for frame in frames if frame["exptime"] < float(floor)]
    keep = [frame for frame in frames if frame not in short]
    if short and len(keep) >= int(min_frames):
        return keep, short
    return list(frames), []


def apply_rate_floor(frames, frac, min_frames):
    """Drop frames under `frac` of the median positive rate. Yields the same way."""
    if frac is None or len(frames) < 2:
        return list(frames), []
    rates = [frame["rate"] for frame in frames
             if frame.get("rate") is not None and frame["rate"] > 0]
    if len(rates) < 2:
        return list(frames), []
    ordered = sorted(rates)
    mid = len(ordered) // 2
    med = ordered[mid] if len(ordered) % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])
    bad = [frame for frame in frames
           if frame.get("rate") is not None and frame["rate"] < float(frac) * med]
    keep = [frame for frame in frames if frame not in bad]
    if bad and len(keep) >= int(min_frames):
        return keep, bad
    return list(frames), []


def split_visits(frames, max_dam, pointing_tol=POINTING_TOL_ARCSEC):
    """One star's frames, in time order, split by airmass spread and pointing."""
    ordered = sorted(frames, key=lambda frame: frame["mjd"])
    if not ordered:
        return []
    groups, cur = [], [ordered[0]]
    for this in ordered[1:]:
        airmass = [frame["airmass"] for frame in cur] + [this["airmass"]]
        moved = max(abs(this["ra"] - cur[0]["ra"]), abs(this["dec"] - cur[0]["dec"]))
        if moved > pointing_tol or max(airmass) - min(airmass) > max_dam:
            groups.append(cur)
            cur = []
        cur.append(this)
    groups.append(cur)
    return groups


def group_frames(frames, max_dam, rate_frac, min_frames, floor):
    """Exposure floor, then per star a count-rate floor, then visit splits.

    Returns groups, frames dropped for exposure time, and frames dropped for
    count rate. A group with fewer than min_frames has refused=True.
    """
    kept, exptime_dropped = apply_exptime_floor(frames, floor, min_frames)
    stars = sorted({frame["star"] for frame in kept})
    groups, rate_dropped = [], []
    for star in stars:
        star_frames = [frame for frame in kept if frame["star"] == star]
        star_frames, dropped = apply_rate_floor(star_frames, rate_frac, min_frames)
        rate_dropped.extend(dropped)
        for visit in split_visits(star_frames, max_dam):
            groups.append({"star": star, "frames": visit,
                           "refused": len(visit) < int(min_frames)})
    return groups, exptime_dropped, rate_dropped


def cube_tag(star, night, letter, index, n_for_star):
    tag = f"{star}_{night}_{letter}"
    if n_for_star > 1:
        tag += f"_v{index}"
    return tag


def render_coadd3d(tag, star, night, letter, visit, n_visits, frames):
    airmass = [frame["airmass"] for frame in frames]
    lines = [
        f"# {star}, {night} setup {letter}, visit {visit} of {n_visits}; "
        f"airmass {min(airmass):.2f}-{max(airmass):.2f}",
        "[rdx]",
        "    spectrograph = keck_kcrm",
        "",
        "[reduce]",
        "    [[cube]]",
        "        combine = True",
        f"        output_filename = {tag}.fits",
        "        save_whitelight = True",
        "",
        "spec2d read",
        "filename",
    ]
    for frame in frames:
        name = frame["spec2d"]
        lines.append(name if name.startswith("Science/") else f"Science/{os.path.basename(name)}")
    lines.append("spec2d end")
    return "\n".join(lines) + "\n"


def region_bounds(region):
    left, right = region.split(",")
    return float(left.lstrip(":") or 0), float(right.rstrip(":") or 100)


def cluster_regions(items, tol=REGION_TOL_PCT):
    """Merge (key, region, frames) whose sky windows agree to within `tol` percent."""
    clusters = []
    for key, region, frames in items:
        bounds = region_bounds(region)
        for cluster in clusters:
            other = region_bounds(cluster["region"])
            if max(abs(bounds[0] - other[0]), abs(bounds[1] - other[1])) <= tol:
                cluster["frames"].extend(frames)
                cluster["keys"].append(key)
                break
        else:
            clusters.append({"region": region, "frames": list(frames), "keys": [key]})
    return clusters


def write_regions(path, clusters):
    lines = ["# sky windows measured by regions.py. science.py reads this file."]
    for index, cluster in enumerate(clusters, start=1):
        lines.append(f"cluster {index} {cluster['region']}")
        for frame in cluster["frames"]:
            lines.append(f"  {frame}")
    text = "\n".join(lines) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return text


def read_regions(path):
    clusters = []
    for line in open(path):
        if line.startswith("cluster "):
            _, _, region = line.split(None, 2)
            clusters.append({"region": region.strip(), "frames": []})
        elif line.startswith("  ") and clusters:
            clusters[-1]["frames"].append(line.strip())
    if not clusters or any(not cluster["frames"] for cluster in clusters):
        raise StepError(f"{path} has no sky windows. Run regions.py again.")
    return clusters


def retype_thar(text, lamp_of):
    """Comment FeAr arc rows and mark ThAr rows arc,tilt.

    `lamp_of(filename)` returns 'FeAr', 'ThAr', or None. A file that already
    contains arc,tilt is left untouched: a second pass would find no tilt rows
    and could comment the only remaining arc.
    """
    if "arc,tilt" in text:
        return text, "already retyped"
    out, moved, dropped = [], 0, 0
    for line in text.split("\n"):
        filename = line.split("|")[0].strip().lstrip("#").strip()
        if not filename.endswith(".fits") or line.lstrip().startswith("#"):
            out.append(line)
            continue
        lamp = lamp_of(filename)
        if lamp == "FeAr":
            out.append("# FeAr, unused: this grating is calibrated on ThAr")
            out.append("# " + line)
            dropped += 1
        elif lamp == "ThAr":
            new = re.sub(
                r"\|(\s*)tilt(\s*)\|",
                lambda m: "|" + " " * max(1, len(m.group(1)) - 5)
                + "arc,tilt" + m.group(2) + "|",
                line, count=1)
            if new == line:
                out.append(line)
            else:
                out.append(new)
                moved += 1
        else:
            out.append(line)
    if not moved:
        raise StepError("no ThAr frames were retyped to arc,tilt. "
                        "This grating's line list needs the ThAr lamp as the arc.")
    return "\n".join(out), f"{moved} ThAr frames -> arc,tilt; {dropped} FeAr frames commented out"


def activate(text, keep):
    """Comment science rows whose filename is not in `keep`. Calibration rows stay."""
    out = []
    for line in text.split("\n"):
        bare = line[1:].strip() if line.lstrip().startswith("#") else line
        # A leading '# ' on a data row. Lines that are prose comments stay.
        leading = line.lstrip()
        commented = leading.startswith("#")
        probe = leading[1:].strip() if commented else leading
        parts = [part.strip() for part in probe.split("|")]
        filename = parts[0] if parts else ""
        is_science = (filename.endswith(".fits") and len(parts) > 1
                      and parts[1] == "science")
        if not is_science:
            out.append(line)
        elif filename in keep:
            out.append(probe if commented else line)
        else:
            out.append(line if commented else "# " + line)
    return "\n".join(out)


def render_sens(blue, red, polyorder):
    blue, red = int(blue), int(red)
    if red < 1:
        raise StepError(f"red trim is {red}. PypeIt does trim_gpm[blue:-red], so 0 masks the whole spectrum.")
    return (
        "# IR sensfunc. trim is the counts cut and the full-slice-coverage cut.\n"
        "[sensfunc]\n"
        f"    algorithm = {SENS_ALGORITHM}\n"
        f"    extr = {SENS_EXTR}\n"
        f"    trim_std_pixs = {blue}, {red}\n"
        f"    polyorder = {int(polyorder)}\n"
    )


def sens_names(tag, polyorder):
    order = int(polyorder)
    return (f"{tag}_IR_p{order}cov.sens", f"{tag}_sens_IR_p{order}cov.fits")
